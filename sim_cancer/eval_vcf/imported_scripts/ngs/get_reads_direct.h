/* written by Jonas Behr, Regina Bohnert and Gunnar Raetsch, FML Tuebingen, Germany, 2010 */

#ifndef __GET_READS_DIRECT_H__
#define __GET_READS_DIRECT_H__

#include <vector>
  using std::vector;
#include "read.h"
#include <stdlib.h>
#include "sam.h"
#include "bam.h"
#include <string>
  using std::string;

//static int g_flag_on = 0, g_flag_off = 0;
static int left_flag_mask = strtol((char*) "0x40", 0, 0);
static int right_flag_mask = strtol((char*) "0x80", 0, 0);
static int reverse_flag_mask = strtol((char*) "0x10", 0, 0);

static int subsample = 1000; 


void parse_cigar(const bam1_t* b, CRead* read)
{
	read->left = (b->core.flag & left_flag_mask) >0;
	read->right = (b->core.flag & right_flag_mask) >0;
	read->reverse = (b->core.flag & reverse_flag_mask) >0;

	read->start_pos = b->core.pos+1;
	read->set_strand('0');
	int len = strlen(bam1_qname(b));
	read->read_id = new char[len+1];
	//sprintf(read->read_id, "%s\0", bam1_qname(b));
	sprintf(read->read_id, "%s", bam1_qname(b));

	for (int k = 0; k < b->core.n_cigar; ++k) 
	{
		int op = bam1_cigar(b)[k] & BAM_CIGAR_MASK; // operation
		int l = bam1_cigar(b)[k] >> BAM_CIGAR_SHIFT; // length
		//int op = bam_cigar_op(bam1_cigar(b)[k]); // operation
		//int l = bam_cigar_oplen(bam1_cigar(b)[k]); // length
		//fprintf(stdout, "op:%d l:%d\n", op, l);
		if (op == BAM_CMATCH) 
		{
			if (k==0)
			{
				read->block_lengths.push_back(l);
				read->block_starts.push_back(0);
			}
			else
			{
				int op_prev = bam1_cigar(b)[k-1] & BAM_CIGAR_MASK; 
				int l_prev = bam1_cigar(b)[k-1] >> BAM_CIGAR_SHIFT;
				if (op_prev==BAM_CREF_SKIP)// intron before
				{
					if (read->block_lengths.size()>=1)
					{
						int last_block_start = (*(read->block_starts.end()-1));
						int intron_start = last_block_start+(*(read->block_lengths.end()-1));
						read->block_lengths.push_back(l);
						read->block_starts.push_back(intron_start+l_prev);
					}
					else
					{
						// start of first block was not a match
						read->block_lengths.push_back(l);
						read->block_starts.push_back(0);
					}
				}
				else
				{
					if (read->block_lengths.size()>=1)
						(*(read->block_lengths.end()-1))+=l;
					else
					{
						read->block_lengths.push_back(l);
						read->block_starts.push_back(0);
					}
				}
			}
		}
		else if (op == BAM_CDEL) 
		{
			if (k>0 && read->block_lengths.size()>=1)
				(*(read->block_lengths.end()-1))+=l;
		} 
		else if (op == BAM_CREF_SKIP)//intron
		{}
		else if (op == BAM_CINS)
		{
			//if (k>0 && read->block_lengths.size()>=1)
			//	(*(read->block_lengths.end()-1))+=l;
		}
		else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)
		{
			read->is_clipped = true;
		}
	}
	// parse auxiliary data
	uint8_t* val = bam_aux_get(b, "XS");
	if (val)
	{
		char strand = bam_aux2A(val);
		if (strand != '+' && strand != '-')
			strand = '.'; 
		 read->set_strand(strand); 
	}
	val = bam_aux_get(b, "NM");
	if (val)
	{
		read->mismatches = bam_aux2i(val);
	}
	val = bam_aux_get(b, "H0");
	if (val)
	{
		read->matches = bam_aux2i(val);
	}
	val = bam_aux_get(b, "HI");
	if (val)
	{
		read->multiple_alignment_index = bam_aux2i(val);
	}


    //uint8_t* s = bam1_aux(b);
	//uint8_t* end = b->data + b->data_len; 
	//while (s < end) 
	//{
	//	 uint8_t type, key[2];
	//	 key[0] = s[0]; key[1] = s[1];
	//	 s += 2; type = *s; ++s;
	//	 //fprintf(stdout, "\n%c%c:%c\n", key[0], key[1], type);
	//	 if (type == 'A')
	//	 {
	//		if ( key[0] =='X' && key[1] == 'S')
	//		{
	//			read->set_strand((char) *s);
	//		}
	//	 	++s;
	//	 }
	//	else if (type == 'C')
	//	{ 
	//		if ( key[0] =='H' && key[1] == '0')
	//		{
	//			uint8_t matches = *s;
	//			read->matches = (int) matches;
	//		}
	//		if ( key[0] =='N' && key[1] == 'M')
	//		{
	//			uint8_t mismatches = *s;
	//			read->mismatches = (int) mismatches;
	//		}
	//		if ( key[0] =='H' && key[1] == 'I')
	//		{
	//			uint8_t mai = *s;
	//			read->multiple_alignment_index = (int) mai;
	//		}

	//		++s;
	//	}
	//	else if (type == 'c') { ++s; }
	//	else if (type == 'S') { s += 2;	}
	//	else if (type == 's') { s += 2;	}
	//	else if (type == 'I') { s += 4; }
	//	else if (type == 'i') { s += 4; }
	//	else if (type == 'f') { s += 4;	}
	//	else if (type == 'd') { s += 8;	}
	//	else if (type == 'Z') { ++s; }
	//	else if (type == 'H') { ++s; }
	//}

	//if (read->strand[0]=='0' && strand_from_flag)
	//{
	//	if ((b->core.flag & reverse_flag_mask) >0)
	//		read->set_strand('-');
	//	else
	//		read->set_strand('+');
	//}
}

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	int beg, end;
	samfile_t *in;
} tmpstruct_t;

typedef struct {
    uint64_t u, v;
} pair64_t;

typedef struct {

	unsigned long num_reads; 
	vector<CRead>* reads; 

} read_buf; 

typedef struct {

	unsigned long num_reads; 
	vector<const bam1_t*>* reads; 

} read_buf_bam1_t; 

static inline int is_overlap(uint32_t beg, uint32_t end, const bam1_t *b)
{
    uint32_t rbeg = b->core.pos;
    uint32_t rend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + 1;
    return (rend > beg && rbeg < end);
}

// callback function for bam_fetch
static int fetch_func(const bam1_t* b, void* data)
{
	read_buf* rbuf = (read_buf*) data; 

	if (rbuf->reads->size()<rbuf->num_reads+1)
	{
		int num_new;
		if (rbuf->num_reads<1e6)
			num_new = 10000+rbuf->num_reads*2;
		else if (rbuf->num_reads<1e7)
			num_new = rbuf->num_reads*1.5;
		else 
			num_new = rbuf->num_reads*1.2;

		rbuf->reads->resize(num_new);
	}
	CRead* read = &(rbuf->reads->at(rbuf->num_reads));
	parse_cigar(b, read);
	rbuf->num_reads++; 

	return 0;
}

static int fetch_bam1_t(const bam1_t* b, void* data)
{
	read_buf_bam1_t* rbuf = (read_buf_bam1_t*) data; 

	if (rbuf->reads->size()<rbuf->num_reads+1)
	{
		int num_new;
		if (rbuf->num_reads<1e6)
			num_new = 10000+rbuf->num_reads*2;
		else if (rbuf->num_reads<1e7)
			num_new = rbuf->num_reads*1.5;
		else 
			num_new = rbuf->num_reads*1.2;

		rbuf->reads->resize(num_new);
	}
	rbuf->reads->at(rbuf->num_reads) = bam_dup1(b); 
	rbuf->num_reads++; 

	return 0;
}

pair64_t * get_chunk_coordinates(const bam_index_t *idx, int tid, int beg, int end, int* cnt_off);

int bam_fetch_reads(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_header_t* header, vector<CRead*>* reads, char strand);

// callback for bam_plbuf_init()
//static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
//{
//	//tmpstruct_t *tmp = (tmpstruct_t*)data;
//	//if ((int)pos >= tmp->beg && (int)pos < tmp->end)
//	//	printf("%s\t%d\t%d\n", tmp->in->header->target_name[tid], pos + 1, n);
//	return 0;
//}

#ifdef __cplusplus
}
#endif

//static int collapse = 0;

//static bool strand_from_flag = false;

int get_reads_from_bam(char* filename, char* region, vector<CRead>* reads, char strand, int lsubsample)
{
	subsample = lsubsample;

	if (strand=='+')
	{}// parameter currently not used, keep it for compatibility reasons
	//set_strand(strand);

	srand (time(NULL));
	//srand (1234);
	tmpstruct_t tmp;
	tmp.in = samopen(filename, "rb", 0);
	if (tmp.in == 0) {
		fprintf(stderr, "Fail to open BAM file %s\n", filename);
		return 1;
	}
	int ref;
	bam_index_t *idx;
	//bam_plbuf_t *buf;
	idx = bam_index_load(filename); // load BAM index
	if (idx == 0) {
		fprintf(stderr, "BAM indexing file is not available.\n");
		samclose(tmp.in);
		return 1;
	}
	bam_parse_region(tmp.in->header, region, &ref,
	                 &tmp.beg, &tmp.end); // parse the region
	if (ref < 0) {
		fprintf(stderr, "Invalid region %s\n", region);
		bam_index_destroy(idx);
		samclose(tmp.in);
		return 1;
	}
	read_buf rbuf; 
	rbuf.reads = reads; 
	rbuf.num_reads = reads->size(); 
	bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, &rbuf, fetch_func);
	reads->resize(rbuf.num_reads);

	bam_index_destroy(idx);
	samclose(tmp.in);
	return 0;
}

int get_reads_from_bam(char* filename, char* region, vector<const bam1_t*>* reads)
{
	tmpstruct_t tmp;
	tmp.in = samopen(filename, "rb", 0);
	if (tmp.in == 0) {
		fprintf(stderr, "Fail to open BAM file %s\n", filename);
		return 1;
	}
	int ref;
	bam_index_t *idx;
	idx = bam_index_load(filename); // load BAM index
	if (idx == 0) {
		fprintf(stderr, "BAM indexing file is not available.\n");
		samclose(tmp.in);
		return 1;
	}
	bam_parse_region(tmp.in->header, region, &ref,
	                 &tmp.beg, &tmp.end); // parse the region
	if (ref < 0) {
		fprintf(stderr, "Invalid region %s\n", region);
		bam_index_destroy(idx);
		samclose(tmp.in);
		return 1;
	}
	read_buf_bam1_t rbuf; 
	rbuf.reads = reads; 
	rbuf.num_reads = reads->size(); 
	bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, &rbuf, fetch_bam1_t);
	reads->resize(rbuf.num_reads);

	bam_index_destroy(idx);
	samclose(tmp.in);
	return 0;
}

#endif
