#include <stdio.h>
#include <vector>
	using std::vector;
#include <string.h>
#include <stdlib.h>
#include "assert.h"
#include <math.h>

#include <ctype.h> // toupper


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
		printf("usage: %s <fn_txt> <fn_vcf> [--max-pval <float>]\n", args[0]); 
		return -1;  
	}
	char* fn_txt = args[1]; 
	char* fn_vcf = args[2]; 

	float filt_pval=10.0; 
	bool filt = false; 
	for (int i=3; i<argc; i++)
	{
		if (strcmp(args[i], "--max-pval")==0)
		{
			assert(i+1<argc); i++;
			filt=true; 
			filt_pval=atof(args[i]); 
		}
		else 
		{
			printf("did not recognize argument: %s\n", args[i]); 
			exit(0); 
		}
	}

	FILE* fd_in = fopen(fn_txt, "r"); 
	assert(fd_in); 	
	FILE* fd_out = fopen(fn_vcf, "w"); 
	assert(fd_out); 

	fprintf(fd_out, "##fileformat=VCFv4.0\n"); 
	fprintf(fd_out, "##source=deepSNV\n"); 
	fprintf(fd_out, "##reference=$fn_genome\n"); 
	fprintf(fd_out, "##INFO=<ID=NV,Number=1,Type=Integer,Description=\"Number of reads supporting variant in normal\">\n"); 
	fprintf(fd_out, "##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Coverage in normal\">\n"); 
	fprintf(fd_out, "##INFO=<ID=TV,Number=1,Type=Integer,Description=\"Number of reads supporting variant in tumor\">\n"); 
	fprintf(fd_out, "##INFO=<ID=CT,Number=1,Type=Integer,Description=\"Coverage in tumor\">\n"); 
	fprintf(fd_out, "##INFO=<ID=VF,Number=.,Type=Float,Description=\"Relative frequency of the SNV\">\n"); 
	fprintf(fd_out, "##INFO=<ID=SF,Number=.,Type=Float,Description=\"Estimated variance of the frequency\">\n"); 
	fprintf(fd_out, "##INFO=<ID=RP,Number=.,Type=Float,Description=\"Raw p-value\">\n"); 
	fprintf(fd_out, "##INFO=<ID=AP,Number=.,Type=Float,Description=\"Adjusted p-value\">\n"); 
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

		if (filt && atof(fields[4])>filt_pval)
			continue; 

		char* chr=fields[0];
		char* pos=fields[1];
		char* ref=fields[2];
		char* var=fields[3];
		char* raw_pVal=fields[4];
		char* freq_var=fields[5] ;
		char* sigma2_freq_var=fields[6];
		char* n_tst_fw=fields[7];
		char* cov_tst_fw=fields[8];
		char* n_tst_bw=fields[9];
		char* cov_tst_bw=fields[10];
		char* n_ctrl_fw=fields[11];
		char* cov_ctrl_fw=fields[12];
		char* n_ctrl_bw=fields[13];
		char* cov_ctrl_bw=fields[14];
		char* adjusted_pVal=fields[15];

		int n_tst=atoi(n_tst_fw) + atoi(n_tst_bw); 
		int cov_tst=atoi(cov_tst_fw) + atoi(cov_tst_bw); 

		int n_ctrl=atoi(n_ctrl_fw) + atoi(n_ctrl_bw); 
		int cov_ctrl=atoi(cov_ctrl_fw) + atoi(cov_ctrl_bw);

		float qual = -10 * log(atof(raw_pVal))/log(2); 
		float VAF_normal = (float(n_ctrl)/cov_ctrl);
		float VAF_tumor = (float(n_tst)/cov_tst);
		float VAF_diff = VAF_tumor - VAF_normal;

		if( VAF_diff >= 0 ){ // then it's a somatic mutation - otherwise it is a LOH event
			fprintf(fd_out, "%s\t%s\t.\t%s\t%s\t%.2f\t.\tNV=%i;CN=%i;TV=%i;CT=%i;VF=%s;SF=%s;RP=%s;AP=%s\t.\n", chr, pos,   ref, var, qual, n_ctrl,  cov_ctrl, n_tst, cov_tst,         freq_var,sigma2_freq_var, raw_pVal, adjusted_pVal); 
		}
		cnt2++;
	}
	printf("\nread %i lines, wrote %i lines, filt_pval: >%.3f\n", cnt, cnt2, filt_pval); 
}


