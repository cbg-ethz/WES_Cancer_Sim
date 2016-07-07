#include <string>
	using std::string;
#include <vector>
	using std::vector;
#include <math.h>
#include <algorithm>
#include <assert.h>
#include "bam.h"
#include "region.h"
#include "get_reads_direct.h"
#include "tools.h"
#include <time.h>
//#include "genome.h"
//#include "gene.h"
//#include "gene_tools.h"
//#include "infer_genes.h"
//#include "Config.h"
#include "string_utils.h"

//includes for samtools
////////////////////////////////////////////////////////////////////////////////////////////////////
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "razf.h"
#include "bgzf.h"
#include "razf.h"
////////////////////////////////////////////////////////////////////////////////////////////////////

struct Config{
	
	bool cache_all; 
	bool no_mm; 
	bamFile fd_in; 
	bamFile fd_out; 
	vector<bam1_t*> all_reads; 
	vector<bam1_t*> reads_out; 
}; 

int get_num_mismatches(bam1_t* b)
{
	//uint8_t* aux = bam_aux_get(reads[i], "XS");
	//char strand = '.'; 
	//int mismatches = 100; 
	//int multi_mappers = 0; 
	//if (aux)
	//	strand = bam_aux2A(aux);
	
	int mismatches = 1000; 
	uint8_t* aux = bam_aux_get(b, "NM");
	if (aux)
		mismatches = bam_aux2i(aux);
	//else
	//	printf("no aux NM: %s\n", bam1_qname(b));
	
	return mismatches; 
	//aux = bam_aux_get(reads[i], "NH");
	//if (aux)
	//	multi_mappers = bam_aux2i(aux);
	
	//aux = bam_aux_get(reads[i], "ZF");
	//if (aux)
	//	mass = bam_aux2i(aux);
}

bool sort_by_ID(bam1_t* b1, bam1_t* b2)
{
	char* id1 = bam1_qname(b1);
	char* id2 = bam1_qname(b2);
	
	return std::lexicographical_compare(id1,id1+strlen(id1),id2,id2+strlen(id2));
}

bool sort_by_pos(bam1_t* b1, bam1_t* b2)
{
	if (b1->core.tid<b2->core.tid)
		return true; 
	if (b1->core.tid>b2->core.tid)
		return false; 

	return (b1->core.pos<b2->core.pos); 
}

void report(Config* c, bam1_t* left, bam1_t* right)
{
	if (c->cache_all)
	{
		if (left)
			c->reads_out.push_back(left); 
		if (right)
			c->reads_out.push_back(right); 
	}
	else
	{
		if (left)
			bam_write1(c->fd_out, left); 
		if (right)
			bam_write1(c->fd_out, right); 
	}
}
void parse_mapping(vector<int*>* genome_map, bam_header_t* header, char* fn_map) 
{
	*genome_map = vector<int*>(header->n_targets, NULL); 

	FILE* fd = fopen(fn_map, "r"); 
	assert(fd); 
	char line[1000]; 
	char chr_prev[100]; 
		sprintf(chr_prev, "xxxxx"); 
	int tid = -1; 
	int cur_tu = 0; 
	int cur_ref = 0; 
	int pos1_prev = 0; 
	int pos2_prev = 0; 
	while (fgets(line, 1000, fd))
	{
		vector<char*> fields = my_str_tok(line, "\t\n"); 
		assert(fields.size()==3); 
		char* chr = fields[0]; 
		int pos1 = atoi(fields[1]); 
		int pos2 = atoi(fields[2]); 

		if (strcmp(chr, chr_prev)!=0)
		{
			// start with new chromosome
			printf("\rstart with chr %s (prev: %s)                                  ", chr, chr_prev); 
			tid = -1; 
			for (int i=0; i<header->n_targets; i++)
			{
				if (strcmp(chr, header->target_name[i])==0)
				{
					sprintf(chr_prev, "%s", chr); 
					tid = i; 
					break; 
				}
			}
			assert(tid>-1);
			assert(genome_map->at(tid)==NULL); 
			genome_map->at(tid) = new int[header->target_len[tid]+1000000]; // the tumor chromosome might be a bit larger
			memset(genome_map->at(tid), 0, (header->target_len[tid]+1000000)*sizeof(int)); 
			cur_tu = 0; 
			cur_ref = 0; 
			pos1_prev = 0; 
			pos2_prev = 0; 
		}
		int posdiff = (pos2-pos2_prev) - (pos1-pos1_prev); 
		
		if (posdiff<0)
		{
			// deleted in tumor
			cur_ref -= posdiff; 
		}
		else if (posdiff>0)
		{
			//inserted in tumor
			for (int i=0; i<posdiff; i++)
			{
				//printf("%i %i\n", cur_ref, cur_tu); 
				genome_map->at(tid)[cur_tu] = cur_ref;
				cur_tu++; 
			}
		}

		while (cur_tu<pos2)
		{
			//printf("%i %i\n", cur_ref, cur_tu); 
			genome_map->at(tid)[cur_tu] = cur_ref; 
			cur_ref++; 
			cur_tu++; 
		}
		assert(cur_ref==pos1); 
		pos1_prev = pos1; 
		pos2_prev = pos2; 

		//if (pos1>10322)
		//	exit(0); 
	}

	fclose(fd); 
}

int main(int argc, char* args[])
{

	if (argc<3)
	{
		printf("Usage: %s <fn_bam> <fn_bam_out> [options]\n", args[0]);
		printf("options:\n"); 
		printf("\t-cache-all\tkeep all reads in memory to sort; otherwise input needs to be sorted \n"); 
		printf("\t-no-mm\tdiscard all reads with multiple mapping locations\n"); 
		exit(1); 
	}

	Config c; 
	c.cache_all = false;
	c.no_mm = false;

	char* fn_bam = args[1]; 
	char* fn_bam_out = args[2]; 
	char* fn_map = NULL; 
	for (int i=3; i<argc; i++)
	{
		if (strcmp(args[i], "--cache-all")==0)
		{
			c.cache_all=true; 
		}
		else if (strcmp(args[i], "--no-mm")==0)
		{
			c.no_mm=true; 
		}
		else if (strcmp(args[i], "--map")==0)
		{
			assert(i<argc-1); 
			i++; 
			fn_map = args[i]; 
		}
		else
		{
			printf("did not understand arg: %s\n", args[i]); 
			exit(1);
		}
	}

	assert(strcmp(fn_bam, fn_bam_out)!=0); 
	printf("open bam file: %s\n", fn_bam); 

	c.fd_in = bam_open(fn_bam, "r");
	if (c.fd_in==0)
	{
		fprintf(stderr, "[%s] Could not open bam file: %s", args[0], fn_bam);
		exit(1);
	}
	bam_header_t* header = bam_header_read(c.fd_in);
	if (header == 0)
	{
		fprintf(stderr, "[%s] Invalid BAM header.", args[0]);
		exit(1);
	}
	c.fd_out = bam_open(fn_bam_out, "w");
	if (c.fd_out==0)
	{
		fprintf(stderr, "[%s] Could not open bam file for writing: %s", args[0], fn_bam_out);
		exit(1);
	}

	vector<int*> genome_map; 
	parse_mapping(&genome_map, header, fn_map); 

	printf("number of chromosomes in header:%i\n", header->n_targets);

	bam_header_write(c.fd_out, header);

	bam1_t* b = NULL;

	unsigned int cnt=0;
	unsigned int cnt_align=0; 
	unsigned int num_same_chr=0; 
	unsigned int cnt_right=0; 
	unsigned int unmap = 0;
	unsigned long mem = 0; 
	int num_read=1;


	if (c.cache_all)
	{
		time_t t1 = time(NULL); 
		int cnt = 0; 
		while(num_read>0)
		{
			bam1_t* b = (bam1_t*) calloc(1, sizeof(bam1_t));
			assert(b);
			num_read = bam_read1(c.fd_in, b); 

			if (num_read>0 && (b->core.flag&BAM_FUNMAP)==0)
			{
				c.all_reads.push_back(b); 
				mem+= b->data_len*sizeof(uint8_t); 
				mem+= sizeof(bam1_t); 
				mem+= sizeof(void*); 
			}
			if (++cnt%1000000==0)
			{
				printf("read %i alignments (mem: %luMb); %.0fsec)\n", cnt, mem/(1024*1024), difftime(time(NULL), t1)); 
			}
		}
		t1 = time(NULL);
		printf("sort by ID...\n"); 
		sort(c.all_reads.begin(), c.all_reads.end(), sort_by_ID); 
		printf("took %.0f seconds\n", difftime(time(NULL), t1));
		num_read = 1;
		printf("read %i alignments (mem: %luMb); %.0fsec)\n", cnt, mem/(1024*1024), difftime(time(NULL), t1)); 
	}
	else
	{
		printf("%s: assume reads are sorted by ID\n", args[0]); 
	}

	// iterate over all raeds
	////////////////////////////////////////////////////////////////////////////////
	while(num_read>0)
	{

		// collect reads with the same read id in vectors reads_left and reads_right
		////////////////////////////////////////////////////////////////////////////////
		vector<bam1_t*> reads_left;  
		vector<bam1_t*> reads_right;  
		char* read_id = NULL; 

		if (b!=NULL)
		{
			if (b->core.flag&BAM_FREAD1)
				reads_left.push_back(b); 
			else
				reads_right.push_back(b); 
			read_id = bam1_qname(b);
		}

		for (;;)
		{
			bam1_t* b2; 
			if (c.cache_all)
			{
				if (cnt_align>= c.all_reads.size()) 
				{
					num_read=0; 
					break; 
				}
				b2 = c.all_reads[cnt_align]; 
			}
			else
			{
				b2 = (bam1_t*) calloc(1, sizeof(bam1_t));
				if ((num_read = bam_read1(c.fd_in, b2))<1) break; 
			}
			bam1_core_t* core2 = &b2->core;

			cnt_align++;
	
			if (core2->flag&BAM_FUNMAP)
			{
				// this read is unaligned
				if (!c.cache_all)
				{
					free(b2->data); 
					free(b2); 
				}
				continue;
				unmap++; 
			}

			if (!read_id)
			{
				read_id = bam1_qname(b2);
				if (b2->core.flag&BAM_FREAD1)
					reads_left.push_back(b2); 
				else
					reads_right.push_back(b2);
				continue;
			}
			
			char* read_id2 = bam1_qname(b2);
			if (strcmp(read_id, read_id2)==0)
			{
				if (b2->core.flag&BAM_FREAD1)
					reads_left.push_back(b2); 
				else
					reads_right.push_back(b2);
				continue;
			}
			
			b = b2; 
			break; 
		}

		if (reads_left.size()==0 && reads_right.size()==0)
			continue; 

		char read_id_cp[1000]; 
		sprintf(read_id_cp, "%s", read_id); 
		vector<char*> fields = my_str_tok(read_id_cp, ";:"); 
		assert(fields.size()==3); 
		vector<char*> fields2 = my_str_tok(fields[1], "_"); 
		if (reads_left.size()==0)
		{
			fields2 = my_str_tok(fields[2], "_"); 
			reads_left = reads_right;
			reads_right.clear(); 
		}
		assert(fields2.size()>=2); 
		char* true_chr = fields2[0]; 
		int true_pos = atoi(fields2[1]); 

		cnt++; 
		//printf("chr: %s pos: %i\n", true_chr, true_pos); 
		bool same_chr=false; 
		int left_idx = -1; 
		int right_idx = -1; 
		for (unsigned int i=0; i<reads_left.size(); i++)
		{
			if (strcmp(true_chr, header->target_name[reads_left[i]->core.tid])==0)
			{
				int true_pos_ref_coord; 
				if (genome_map[reads_left[i]->core.tid])
					true_pos_ref_coord = genome_map[reads_left[i]->core.tid][true_pos]; 
				else 
					true_pos_ref_coord = true_pos; 

				if (abs(reads_left[i]->core.pos - true_pos_ref_coord)<30)
				{
					same_chr = true; 
					left_idx = i; 
					break; 
				}
				//else if (abs(reads_left[i]->core.pos - true_pos_ref_coord)<130)
				//{
				//	printf("align_pos: %i true_pos:%i diff %i\n", reads_left[i]->core.pos, true_pos_ref_coord, reads_left[i]->core.pos - true_pos_ref_coord); 
				//	assert(false); 
				//}
			}

			//if (same_chr)
			//	printf("same: %s %s\n", true_chr, header->target_name[reads_left[i]->core.tid]); 
			//else
			//	printf("diff: %s %s\n", true_chr, header->target_name[reads_left[i]->core.tid]); 
		}
		num_same_chr += same_chr; 
		if (cnt%1000000==0)
			printf("num_reads:%i num same: %i (%.6f%%)\n", cnt, num_same_chr, ((float) num_same_chr)/cnt); 

		if (!same_chr)
			continue; 
		
		// find the matching right read 
		for (unsigned int i=0; i<reads_right.size(); i++)
		{
			if (reads_right[i]->core.pos==reads_left[left_idx]->core.mpos)
			{
				right_idx = i; 
			}
		}

		if (right_idx>=0)
		{
			report(&c, reads_left[left_idx], reads_right[right_idx]);
			cnt_right++; 
		}
		else
			report(&c, reads_left[left_idx], NULL);
	}
	if (c.cache_all)
	{
		printf("sort %lu output reads\n", c.reads_out.size()); 
		sort(c.reads_out.begin(), c.reads_out.end(), sort_by_pos); 

		for (unsigned int i=0; i<c.reads_out.size(); i++)
		{
			bam_write1(c.fd_out, c.reads_out[i]);
		}
		for (unsigned int i=0; i<c.all_reads.size(); i++)
			free(c.all_reads[i]);
	}

	bam_close(c.fd_in);
	bam_close(c.fd_out);

	printf("found %i reads, %u alignments; kept:%i (right:%i) reads (%.6f%%)", cnt, cnt_align, num_same_chr, cnt_right, ((float) num_same_chr)/cnt); 
}
	
