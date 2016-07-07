#include <stdio.h>
#include <vector>
	using std::vector;
#include <string.h>
#include <stdlib.h>
#include "assert.h"


//std::vector<char*> my_str_tok(char* line, const char* sep)
void my_str_tok(char* line, const char* sep, vector<char*>* ret)
{
//	std::vector<char*> ret;
	char* tok = strtok(line, sep); 
	while (tok)
	{
		ret->push_back(tok); 
		tok = strtok(NULL, sep);
	}
	//return ret;
}

int main(int argc, char** args)
{
	if (argc<2)
	{
		printf("usage: %s <fn_txt> <fn_vcf> [--min-prob <float>]\n", args[0]); 
		return -1;  
	}
	char* fn_txt = args[1]; 
	char* fn_vcf = args[2]; 

	float filt_prob=-1;
	bool filt = false; 
	if (argc>=4)
	{
		if (strcmp(args[3], "--min-prob")==0)
		{
			filt=true; 
			filt_prob=atof(args[4]); 
		}
	}

	FILE* fd_in = fopen(fn_txt, "r"); 
	assert(fd_in); 	
	FILE* fd_out = fopen(fn_vcf, "w"); 
	assert(fd_out); 

	fprintf(fd_out, "##fileformat=VCFv4.0\n"); 
	fprintf(fd_out, "##source=JointSNVMix\n"); 
	fprintf(fd_out, "##reference=$fn_genome\n"); 
	fprintf(fd_out, "##INFO=<ID=NR,Number=1,Type=Integer,Description=\"Number of reads supporting reference in normal\">\n"); 
	fprintf(fd_out, "##INFO=<ID=NV,Number=1,Type=Integer,Description=\"Number of reads supporting variant in normal\">\n");
	fprintf(fd_out, "##INFO=<ID=TR,Number=1,Type=Integer,Description=\"Number of reads supporting reference in tumor\">\n");
	fprintf(fd_out, "##INFO=<ID=TV,Number=1,Type=Integer,Description=\"Number of reads supporting variant in tumor\">\n"); 
	fprintf(fd_out, "##INFO=<ID=RPS,Number=.,Type=Float,Description=\"Raw probability of somatic mutation P=p_AA_AB+p_AA_BB\">\n");
	fprintf(fd_out, "##INFO=<ID=PS,Number=.,Type=Float,Description=\"Post processed probability of somatic mutation\">\n");
	fprintf(fd_out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"); 

	char line[10000];
	fgets(line, 10000, fd_in); // skip header

	int cnt=0; 
	int cnt2=0; 
	while (fgets(line, 10000, fd_in)) 
	{
		if (++cnt%100000==0)
			printf("\rline: %i", cnt); 

		vector<char*> fields; 
		my_str_tok(line, "\t\n\"", &fields);

		if (filt && atof(fields[17])<filt_prob)
			continue; 

		char* chr=fields[0];
		char* pos=fields[1];
		char* ref=fields[2];
		char* var=fields[3];
		char* normal_a=fields[4];
		char* normal_b=fields[5] ;
		char* tumor_a=fields[6];
		char* tumor_b=fields[7];

		char* p_AA_AB=fields[9];
		char* p_AA_BB=fields[10];
		char* post_processed_p_somatic=fields[17];
		float post_processed_p_somatic_float=atof(post_processed_p_somatic);

		float PS=atof(p_AA_AB) + atof(p_AA_BB); 

		fprintf(fd_out, "%s\t%s\t.\t%s\t%s\t%f\t.\tNR=%s;NV=%s;TR=%s;TV=%s;RPS=%f;PS=%f\t.\n", chr, pos,   ref, var, post_processed_p_somatic_float,  normal_a,normal_b, tumor_a, tumor_b,PS,post_processed_p_somatic_float); 

		cnt2++;
	}
	printf("\nread %i lines, wrote %i lines, filt_prob: >%.3f\n", cnt, cnt2, filt_prob); 
}


