#ifndef __BAM_TOOLS_H__
#define __BAM_TOOLS_H__

#include <vector>
	using std::vector;
#include "ngs/read.h"
#include "ngs/get_reads_direct.h"

#ifdef _REGIONS_H__
// avoid dependency on regions.h for tools that dont use it
//#include "ngs/regions.h"
void get_regions(vector<Region*>* regions, uint32_t* map, int len, char* chr, char strand, uint16_t shrink)
{
	int map_resolution = 1; 

	bool in_reg = false;
	int start = -1;
	int stop = -1;
	for (int i=0; i<len; i++)
	{
		if (!in_reg && map[i]>0)
		{
			start = i*map_resolution;
			stop = -1;
			in_reg = true;
		}
		else if (in_reg && map[i]==0)
		{
			stop = i*map_resolution;
			if (shrink>1)
			{
				//fprintf(stdout, "start:%i, stop:%i", start, stop);
				// move start
				for (int st=start/map_resolution; st<len && map[st]<shrink; st++)
					start = st*map_resolution;
				// move stop
				for (int st=i; st>start/map_resolution && map[st]<shrink; st--)
					stop = st*map_resolution;
				//fprintf(stdout, "  shrink-> start:%i, stop:%i\n", start, stop);
			}

			if (stop-start>1)
			{
				char* pchr = new char[strlen(chr)+1];
				strcpy(pchr, chr);
				Region* reg = new Region(start, stop, pchr, strand);
				regions->push_back(reg);
			}
			start = -1;
			stop = -1;
			in_reg = false;
		}
	}
}
#endif

void bam_compute_coverage(char* fn_bam, char* chr, int start, int stop, vector<uint32_t>* coverage, bool include_introns) 
{
	int num_pos = stop-start+1;
	if ((int) coverage->size()!=num_pos)
	{
		assert(coverage->size()==0); 
		*coverage = vector<uint32_t>(num_pos, 0);
	}

	vector<CRead> reads; 

	char reg_str[1000]; 
	sprintf(reg_str, "%s:%i-%i", chr, start, stop); 
	
	char strand='.'; 
	int subsample = 1000; // no subsampling
	get_reads_from_bam(fn_bam, reg_str, &reads, strand, subsample); 

	//memset(coverage, 0, num_pos*sizeof(uint32_t));
	//for (int i=0; i<num_pos; i++)
	//	coverage[i] = 0;
	for (uint i=0; i<reads.size(); i++)
	{
		// add read contribution to coverage
		if (include_introns)
			reads[i].get_coverage_all(start, stop, &(coverage->at(0)));
		else
			reads[i].get_coverage(start, stop, &(coverage->at(0)));
	}
}
void bam_compute_coverage(char* fn_bam, char* chr, int start, int stop, vector<uint32_t>* coverage) 
{
	bam_compute_coverage(fn_bam, chr, start, stop, coverage, false); 
}

vector<char*> bam_get_chr_names(char* fn_bam) 
{
	vector<char*> ret; 

	bamFile fd1 = bam_open(fn_bam, "r");
	if (fd1==0)
	{   
		fprintf(stderr, "[%s] Could not open bam file: %s", __func__, fn_bam);
		exit(-1);
	}
	bam_header_t* header = bam_header_read(fd1);


	for (int i=0; i<header->n_targets; i++)
	{
		char* name = new char[strlen(header->target_name[i])+1]; 
		sprintf(name, "%s",  header->target_name[i]); 
		ret.push_back(name); 
	}
	bam_header_destroy(header);
	bam_close(fd1);
	return ret; 
}

uint32_t bam_get_chr_len(char* fn_bam, char* chr) 
{
	bamFile fd1 = bam_open(fn_bam, "r");
	if (fd1==0)
	{   
		fprintf(stderr, "[%s] Could not open bam file: %s", __func__, fn_bam);
		exit(-1);
	}
	bam_header_t* header = bam_header_read(fd1);


	for (int i=0; i<header->n_targets; i++)
	{   
		if (strcmp(chr, header->target_name[i])==0)
		{
			uint32_t len = header->target_len[i]; 
			bam_header_destroy(header);
			bam_close(fd1);
			return len; 
		}
	}
	printf("could not find chromosome: %s in bam header %s\n", chr, fn_bam); 
	printf("chr names are:\n");
	for (int i=0; i<header->n_targets; i++)
	{   
		printf("%s ", header->target_name[i]); 
	}
	printf("\n"); 
	bam_header_destroy(header);
	bam_close(fd1);

	exit(0); 
	return 0; 
}

#endif
