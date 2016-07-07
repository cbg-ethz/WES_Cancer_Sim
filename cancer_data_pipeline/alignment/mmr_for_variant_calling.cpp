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

#define WRITE_MIN_PAIRS

struct Config{
	
	bool cache_all; 
	bool no_mm; 
	bamFile fd_in; 
	bamFile fd_out; 
	vector<bam1_t*> all_reads; 
	vector<bam1_t*> reads_out; 
#ifdef WRITE_MIN_PAIRS
	FILE* fd_min_pairs; 
#endif

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

bool report_single(Config* c, vector<bam1_t*>* reads, int min_mm, int min_cnt)
{

	if (min_cnt==1)
	{
		for (int i=0; i<reads->size(); i++)
		{
			int mm = get_num_mismatches(reads->at(i)); 
			if (mm==min_mm)
			{
				report(c, reads->at(i), NULL);
				return true; 
			}
		}
		assert(false); 
	}

	// choose any of the reads at random
	int idx = rand()%min_cnt; 
	min_cnt=0;
	for (int i=0; i<reads->size(); i++)
	{
		int mm = get_num_mismatches(reads->at(i)); 
		if (mm==min_mm && min_cnt==idx)
		{
			report(c, reads->at(i), NULL);
			return true;
		}
		else if (mm==min_mm)
			min_cnt++;
	}
	return false; 
}

bool  report_no_mm(Config* c, vector<bam1_t*>* reads1, vector<bam1_t*>* reads2)
{
	// report only reads that have no multi mappers
	if (reads1->size()==1 && reads2->size()==1)
		report(c, reads1->at(0), reads2->at(0)); 
	else if (reads1->size()==1 && reads2->size()==0)
		report(c, reads1->at(0), NULL); 
	else if (reads1->size()==0 && reads2->size()==1)
		report(c, NULL, reads2->at(0)); 
	else 
		return false;

	return true; 
}

bool report_min_mm(Config* c, vector<bam1_t*>* reads1, vector<bam1_t*>* reads2, int min_mm1, int min_mm2, bool take_any)
{

	// count the number of pairs that have minimal number of mismatches
	int cnt = 0; 
	for (int i=0; i<reads1->size(); i++)
	{
		int mm1 = get_num_mismatches(reads1->at(i));
		if (mm1!=min_mm1)
		{
			continue; 
		}
		for (int j=0; j<reads2->size(); j++)
		{
			if (reads2->at(j)->core.pos==reads1->at(i)->core.mpos)
			{
				int mm2 = get_num_mismatches(reads2->at(j));
				if (mm2==min_mm2)
				{
					cnt++; 
				}
				break; 
			}
		}
	}
#ifdef WRITE_MIN_PAIRS
	fprintf(c->fd_min_pairs, "%i\n", cnt); 
#endif

	if (cnt==0)
		return false; 

	if (cnt>1 && !take_any)
		return false; 
	
	int idx = rand()%cnt; 

	cnt = 0; 
	for (int i=0; i<reads1->size(); i++)
	{
		int mm1 = get_num_mismatches(reads1->at(i));
		if (mm1!=min_mm1)
		{
			continue; 
		}
		for (int j=0; j<reads2->size(); j++)
		{
			if (reads2->at(j)->core.pos==reads1->at(i)->core.mpos)
			{
				int mm2 = get_num_mismatches(reads2->at(j));
				if (mm2==min_mm2 && cnt==idx)
				{
					report(c, reads1->at(i), reads2->at(j)); 
					return true; 
				}
				else if (mm2==min_mm2)
				{
					cnt++; 
				}
				break; 
			}
		}
	}

	assert(false); 
	return false; 
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

#ifdef WRITE_MIN_PAIRS
	c.fd_min_pairs = fopen("cnt_min_pairs.txt", "w"); 
	assert(c.fd_min_pairs); 
#endif

	char* fn_bam = args[1]; 
	char* fn_bam_out = args[2]; 
	for (int i=3; i<argc; i++)
	{
		if (strcmp(args[i], "-cache-all")==0)
		{
			c.cache_all=true; 
		}
		else if (strcmp(args[i], "-no-mm")==0)
		{
			c.no_mm=true; 
		}
		else
		{
			printf("did not understand arg: %s\n", args[i]); 
			exit(1);
		}
	}

	assert(strcmp(fn_bam, fn_bam_out)!=0); 

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
	if (c.fd_in==0)
	{
		fprintf(stderr, "[%s] Could not open bam file for writing: %s", args[0], fn_bam_out);
		exit(1);
	}

	printf("number of chromosomes in header:%i\n", header->n_targets);

	bam_header_write(c.fd_out, header);

	bam1_t* b = NULL;

	int cnt=0;
	int cnt_align=0; 
	int discard00 = 0;
	int discard01 = 0;
	int discard10 = 0;
	int discard11 = 0;
	int mm00 = 0; 
	int mm01 = 0; 
	int mm10 = 0; 
	int mm11 = 0; 
	int unmap = 0;
	int num_single = 0; 
	int num_pairs = 0; 
	unsigned long mem=0; 
	int len;
	char* chr; 
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
				printf("read %i alignments (mem: %iMb); %.0fsec)\n", cnt, mem/(1024*1024), difftime(time(NULL), t1)); 
			}
		}
		t1 = time(NULL);
		printf("sort by ID...\n"); 
		sort(c.all_reads.begin(), c.all_reads.end(), sort_by_ID); 
		printf("took %.0f seconds\n", difftime(time(NULL), t1));
		num_read = 1;
		printf("read %i alignments (mem: %iMb); %.0fsec)\n", cnt, mem/(1024*1024), difftime(time(NULL), t1)); 
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
				if (cnt_align>=c.all_reads.size()) 
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
		//printf("left: %lu right: %lu cnt:%i num_read: %i\n", reads_left.size(), reads_right.size(), cnt, num_read); 

		// compute stats for the number of mismatches in the alignments
		////////////////////////////////////////////////////////////////////////////////
		int min_mm_left = 10000; 
		int min_cnt_left = 0;
		for (int i=0; i<reads_left.size(); i++)
		{
			int mm = get_num_mismatches(reads_left[i]); 

			if (mm<min_mm_left)
			{
				min_cnt_left = 1;
				min_mm_left = mm; 
			}
			else if (mm==min_mm_left)
				min_cnt_left++; 
		}

		int min_mm_right = 10000; 
		int min_cnt_right = 0;
		for (int i=0; i<reads_right.size(); i++)
		{
			int mm = get_num_mismatches(reads_right[i]); 

			if (mm<min_mm_right)
			{
				min_cnt_right = 1;
				min_mm_right = mm; 
			}
			else if (mm==min_mm_right)
				min_cnt_right++; 
		}

		if (c.no_mm)
		{
			num_pairs++; 
			num_single += report_no_mm(&c, &reads_left, &reads_right); 
			continue; 
		}
		//if (cnt_align>1000000 && cnt_align<1000999)
		//	printf("%s: mm_left: %i cnt: %i  mm_right: %i, cnt:%i\n",read_id, min_mm_left, min_cnt_left, min_mm_right, min_cnt_right);

		// decide on the alignment pair to report
		////////////////////////////////////////////////////////////////////////////////
		// cases: 
		// for x>1 and k>0
		// id |  left       |    right   |                   action
		//    | min_mm  cnt | min_mm cnt |
		// ----------------------------------------------------------------------------------
		// 1  |   0     x   |   0     x  | take if there is one matching pair with both zeros
		// 2  |   0     x   |   k     1  | take according to right
		// 3  |   0     x   |   k     x  | take only if there is only one such pair
		// 4  |   k     1   |   0     x  | like case two
		// 5  |   k     1   |   k     1  | take if there is a match
		// 6  |   k     1   |   k     x  | take according to left
		// 7  |   k     x   |   0     x  | take only if there is only one such pair
		// 8  |   k     x   |   k     1  | take according to right
		// 9  |   k     x   |   k     x  | take only if there is only one such pair

		if (reads_right.size()==0)
		{
			if (min_mm_left>0 && min_cnt_left>1)
				discard10++; 
			else if (min_mm_left==0)
			{
				mm00++; 
				assert(report_single(&c, &reads_left, min_mm_left, min_cnt_left)); 
			}
			else
			{
				mm10++; 
				assert(report_single(&c, &reads_left, min_mm_left, min_cnt_left)); 
			}
		}
		else if (reads_left.size()==0)
		{
			if (min_mm_right>0 && min_cnt_right>1)
				discard01++; 
			else if (min_mm_right==0)
			{
				mm00++; 
				assert(report_single(&c, &reads_right, min_mm_right, min_cnt_right)); 
			}
			else
			{
				mm01++; 
				assert(report_single(&c, &reads_right, min_mm_right, min_cnt_right)); 
			}
		}
		else if (min_mm_left==0 && min_mm_right==0)
		{
			// case 1
			// 1  |   0     x   |   0     x  | take if there is one matching pair with both zeros
			bool take_any = true;
			bool ret = report_min_mm(&c, &reads_left, &reads_right, min_mm_left, min_mm_right, take_any); 
			if (ret)
				mm00++; 
			else
				discard00++; 
		}
		else if (min_mm_left==0 && min_mm_right>0 )
		{
			// case 2 and 3 
			// 2  |   0     x   |   k     1  | take only if there is only one such pair
			// 3  |   0     x   |   k     x  | take only if there is only one such pair
			bool take_any = false;
			bool ret = report_min_mm(&c, &reads_left, &reads_right, min_mm_left, min_mm_right, take_any); 
			if (ret)
				mm01++; 
			else
				discard01++; 
		}
		else if (min_mm_left>0 && min_cnt_left==1 && min_mm_right==0)
		{
			// case 4
			// 4  |   k     1   |   0     x  | take according to left
			bool take_any = false;
			bool ret = report_min_mm(&c, &reads_left, &reads_right, min_mm_left, min_mm_right, take_any); 
			if (ret)
				mm10++; 
			else
				discard10++; 
		}
		else if (min_mm_left>0 && min_cnt_left==1 && min_mm_right>0)
		{
			// case 5 and 6 
			// 5  |   k     1   |   k     1  | take if there is a match
			// 6  |   k     1   |   k     x  | take if there is a match
			bool take_any = false;
			bool ret = report_min_mm(&c, &reads_left, &reads_right, min_mm_left, min_mm_right, take_any); 
			if (ret)
				mm11++; 
			else
				discard11++; 
		}
		else if (min_mm_left>0 && min_cnt_left>1 && min_mm_right==0)
		{
			// case 7
			// 7  |   k     x   |   0     x  | take only if there is only one such pair
			bool take_any = false;
			bool ret = report_min_mm(&c, &reads_left, &reads_right, min_mm_left, min_mm_right, take_any); 
			if (ret)
				mm10++; 
			else
				discard10++; 
		}
		else if (min_mm_left>0 && min_cnt_left>1 && min_mm_right>0)
		{
			// case 8 and 9
			// 8  |   k     x   |   k     1  | take if pair matches
			// 9  |   k     x   |   k     x  | take only if there is only one such pair
			bool take_any = false;
			bool ret = report_min_mm(&c, &reads_left, &reads_right, min_mm_left, min_mm_right, take_any); 
			if (ret)
				mm11++; 
			else
				discard11++; 
		}
		else
		{
			printf("Error, we missed this case\n"); 
			exit(1); 
		}

		// clear vector
		if (!c.cache_all)
		{
			for (int i=0; i<reads_left.size(); i++)
			{
				free(reads_left[i]->data); 
				free(reads_left[i]); 
			}
			for (int i=0; i<reads_right.size(); i++)
			{
				free(reads_right[i]->data); 
				free(reads_right[i]); 
			}

		}
		reads_left.clear(); 
		reads_right.clear(); 

		cnt++;
		if (cnt%1000000==0)
			printf("read %i reads (%i alignment): discarded:%i(%i,%i,%i,%i), no mm:%i only right mm: %i only left mm:%i both mm:%i\n", cnt, cnt_align, discard00+discard01+discard10+discard11, discard00, discard01, discard10, discard11, mm00, mm01, mm10, mm11); 
	}
	printf("read %i reads (%i alignment): discarded:%i(%i,%i,%i,%i), no mm:%i only right mm: %i only left mm:%i both mm:%i\n", cnt, cnt_align, discard00+discard01+discard10+discard11, discard00, discard01, discard10, discard11, mm00, mm01, mm10, mm11); 
	//printf("read %i reads (%i alignment): discarded:%i, no mm:%i only right mm: %i only left mm:%i both mm:%i\n", cnt, cnt_align, discard, mm00, mm01, mm10, mm11); 

	if (c.cache_all)
	{
		printf("sort %lu output reads\n", c.reads_out.size()); 
		sort(c.reads_out.begin(), c.reads_out.end(), sort_by_pos); 

		for (int i=0; i<c.reads_out.size(); i++)
		{
			bam_write1(c.fd_out, c.reads_out[i]);
		}
		for (int i=0; i<c.all_reads.size(); i++)
			free(c.all_reads[i]);
	}
	if (c.no_mm)
		printf("%i out of %i (%.2f%) reads/read pairs with just a single mapping location\n", num_single, num_pairs, ((float)num_single)/num_pairs); 

	bam_close(c.fd_in);
	bam_close(c.fd_out);
	#ifdef WRITE_MIN_PAIRS
		fclose(c.fd_min_pairs); 
	#endif

}
	
