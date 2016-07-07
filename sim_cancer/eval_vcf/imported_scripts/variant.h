
#ifndef __VARIANT_H__
#define __VARIANT_H__

#include <stdio.h>
#include <algorithm> 
	using std::lexicographical_compare; 
#include <string.h>
#include <vector>
	using std::vector; 
#include <assert.h>

// source 
#define BCF 1
#define SDI 2
#define CNV 3
#define RDIFF 5
#define SOM_SNIPER 6
#define PINDEL 7


// location
#define EXONIC 1
#define INTRONIC 2
#define INTERGENIC 3

struct variant
{
	char* chr; 
	int pos;
	int sample;
	int conf_cnt;
	int non_conf_cnt;
	int conf_cnt_normal; 
	int non_conf_cnt_normal; 
	float freq; 
	float qual;
	int source;
	int location;
	bool synonymous;
	char* ref;
	char* alt;
	char ref_AA; 
	char alt_AA; 
	bool splice_site; 
	char* type; 
	char* gene_name; 
	char* info; 
	char* gt_normal; 
	char* gt_cancer; 
};

variant init_var();
variant copy_var(variant v);
bool variant_equal(variant v1, variant v2);
bool variant_equal_seq(variant v1, variant v2);
bool variant_compare(variant v1, variant v2);
bool variant_compare_seq(variant v1, variant v2);
bool variant_qual(variant v1, variant v2);
int print_vcf(FILE* fd, variant v);
void write_vcf(char* vcf_out, vector<variant>* all_var); 

variant init_var()
{
	variant v; 
	v.chr = NULL; 
	v.pos = -1;
	v.sample = -1;
	v.conf_cnt = -1;
	v.non_conf_cnt = -1;
	v.conf_cnt_normal = -1;
	v.non_conf_cnt_normal = -1;
	v.qual = -1.0;
	v.source = 0;
	v.location = 0;
	v.synonymous = true; 
	v.ref = NULL;
	v.alt = NULL;
	v.ref_AA = '.'; 
	v.alt_AA = '.'; 
	v.splice_site=false;
	v.type = NULL; 
	v.gene_name = NULL; 
	v.info = NULL; 
	v.gt_normal = NULL; 
	v.gt_cancer = NULL; 

	return v;
}

char* cp_str(char* str)
{
	if (str)
	{
		char* ret = new char[strlen(str)+1]; 
		sprintf(ret, "%s", str); 
		return ret; 
	}
	return NULL; 
}

variant copy_var(variant v)
{
	variant ret = init_var(); 
	ret.chr = cp_str(v.chr); 
	ret.pos = v.pos;
	ret.sample = v.sample;
	ret.conf_cnt = v.conf_cnt;
	ret.non_conf_cnt = v.non_conf_cnt;
	ret.conf_cnt_normal = v.conf_cnt_normal;
	ret.non_conf_cnt_normal = v.non_conf_cnt_normal;
	ret.qual = v.qual;
	ret.source = v.source;
	ret.location = v.location;
	ret.synonymous = v.synonymous; 
	ret.ref =  cp_str(v.ref);
	ret.alt =  cp_str(v.alt);
	ret.ref_AA = v.ref_AA; 
	ret.alt_AA = v.alt_AA; 
	ret.splice_site=v.splice_site;
	ret.type = cp_str(v.type); 
	ret.gene_name = cp_str(v.gene_name); 
	ret.info = cp_str(v.info); 
	ret.gt_normal = cp_str(v.gt_normal); 
	ret.gt_cancer = cp_str(v.gt_cancer); 
	return ret; 
}

int print_vcf(FILE* fd, variant v)
{
	//#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  H_blood H_after
	if (!v.info)
		fprintf(fd, "%s\t%i\t.\t%s\t%s\t%.3f\t.\t.\t.\n", v.chr, v.pos, v.ref, v.alt, v.qual); 
	else
		fprintf(fd, "%s\t%i\t.\t%s\t%s\t%.3f\t.\t%s\t.\n", v.chr, v.pos, v.ref, v.alt, v.qual, v.info); 
	return 0; 
}

void write_vcf(char* vcf_out, vector<variant>* all_var) 
{
	FILE* fd = fopen(vcf_out, "w"); 
	assert(fd); 
	fprintf(fd, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"); 
	for (unsigned int i=0; i<all_var->size(); i++)
	{
		print_vcf(fd, all_var->at(i)); 
	}
	fclose(fd); 
}

bool variant_equal(variant v1, variant v2)
{
	if (v1.pos!=v2.pos)
		return false; 
	return strcmp(v1.chr, v2.chr)==0;
}
bool variant_compare(variant v1, variant v2)
{
	if (v1.pos<v2.pos)
		return true; 
	if ((v1.pos>v2.pos))
		return false; 
	
	int len1 = strlen(v1.chr); 
	int len2 = strlen(v2.chr); 

	return lexicographical_compare(v1.chr, v1.chr+len1, v2.chr, v2.chr+len2);
}

bool variant_equal_seq(variant v1, variant v2)
{
	if (v1.pos!=v2.pos)
		return false; 
	if (strcmp(v1.chr, v2.chr)!=0)
		return false; 
	//if (strcmp(v1.ref, v2.ref)!=0)
	if (strlen(v1.ref)!=strlen(v2.ref))// some tools define the reference base a the major allel in the normal 
		return false; 
	if (strcmp(v1.alt, v2.alt)!=0)
		return false; 
	return true; 
}

bool variant_compare_seq(variant v1, variant v2)
{
	if (v1.pos<v2.pos)
		return true; 
	if ((v1.pos>v2.pos))
		return false; 
	if (strcmp(v1.chr, v2.chr)!=0)
		return lexicographical_compare(v1.chr, v1.chr+strlen(v1.chr), v2.chr, v2.chr+strlen(v2.chr));
	if (strcmp(v1.ref, v2.ref)!=0)
		return lexicographical_compare(v1.ref, v1.ref+strlen(v1.ref), v2.ref, v2.ref+strlen(v2.ref));

	return lexicographical_compare(v1.alt, v1.alt+strlen(v1.alt), v2.alt, v2.alt+strlen(v2.alt));
}

bool variant_qual(variant v1, variant v2)
{
	return (v1.qual>v2.qual); 
}
#endif
