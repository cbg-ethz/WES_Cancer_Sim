#ifndef __EVALUATE_VCF_H__
#define __EVALUATE_VCF_H__

#include "stats_tests.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <tuple>
#include <cstdlib>


// categories of false positives
#define CLOSE_TO_INDEL 1
#define MISSED_IN_NORMAL 2
#define NOT_IN_TRUE_ALIGN 4


bool pair_comp(pair<float,int> p1, pair<float,int> p2)
{
	return p1.first<p2.first; 
}

char seqi_to_char(uint8_t si)
{
	// the table is: char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";
	//return bam_nt16_rev_table[si]; 
	//char Mybam_nt16_rev_table[] = "=ACMGRSVTWYHKDBN"; 
	const char myTable[] = "=ACMGRSVTWYHKDBN"; // this is from /samtools-0.1.19/bam_import.c
	//return Mybam_nt16_rev_table[si]; 
	return myTable[si];
}

bool read_supports_variant(const bam1_t* b, variant v, bool* supp_ref)
{
	// parse cigar	
	int start_pos = b->core.pos+1; 
	int len_r=0; // position in read sequence
	int len_g=0; // position in genome sequence
	for (int k = 0; k < b->core.n_cigar; ++k) 
	{
		int op = bam_cigar_op(bam1_cigar(b)[k]); // operation
		int l = bam_cigar_oplen(bam1_cigar(b)[k]); // length
		
		//printf("%i %i %i %i %i \n", BAM_CDEL, BAM_CREF_SKIP, BAM_CINS, BAM_CSOFT_CLIP, BAM_CHARD_CLIP); 
		//printf("start_pos:%i len:%i, v.pos:%i\n", start_pos, len_r, v.pos); 

		if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP )
		{
			len_r+=l; 
			//printf("add r: op: %i, len:%i\n", op, l); 
		}
		else if (op == BAM_CDEL || op == BAM_CREF_SKIP ) // deletion in read or intron
		{
			if (start_pos+len_g + l >=v.pos) 
			{
				// variant is in deletion or intron;
				return false; 
			}
			//printf("add g: op: %i, len:%i\n", op, l); 
			len_g+=l; 
		}
		else 
		{
			if (start_pos+len_g + l >=v.pos) 
			{
				len_r += v.pos - start_pos - len_g ; 
				break;
			}
			//printf("add rg: op: %i, len:%i\n", op, l); 
			len_r+=l; 
			len_g+=l; 
		}
	}

	bool supp_alt = true; 
	for (uint8_t j=0; j<strlen(v.alt); j++)
	{
		//std::cerr << "[" << __func__ << "] v.alt = " << v.alt << "\n";
		uint8_t si = bam1_seqi(bam1_seq(b), len_r+j); 
		//std::cerr << "[" << __func__ << "] si = " << static_cast<int>(si) << "\t" << si << "\n";
		char s = seqi_to_char(static_cast<int>(si)); 
		//char fakeS = seqi_to_char(si);
		//std::cerr << "[" << __func__ << "] s = " << s << "\t" << fakeS << "\n";
		//for(uint8_t lol=0; lol < 16; lol++){
		//	std::cerr << "[" << __func__ << "] lol =" << lol << ", and seqi_to_char(lol)="<<  seqi_to_char(lol) << "\n";
		//}
		if (len_r+j>=b->core.l_qseq || v.alt[j]!=s)
			supp_alt = false; 
	}
	*supp_ref = true; 
	for (uint8_t j=0; j<strlen(v.ref); j++)
	{
		uint8_t si = bam1_seqi(bam1_seq(b), len_r+j); 
		//char s = seqi_to_char(si); 
		char s = seqi_to_char(static_cast<int>(si));
		if (len_r+j>=b->core.l_qseq || v.ref[j]!=s)
			*supp_ref = false; 
	}

	//if (false)
	//{
	//	printf("ref: %s, alt: %s string: ", v.ref, v.alt); 
	//	for (int j=-3; j<((int) strlen(v.ref)+3); j++)
	//	{
	//		uint8_t si = bam1_seqi(bam1_seq(b), len_r+j); 
	//		printf("%c", seqi_to_char(si));
	//	}
	//	printf("supp_ref: %i supp_alt: %i\n", supp_ref, supp_alt); 
	//	
	//}
	return supp_alt; 
}

void check_bam_qual(char* fn_bam, variant v, int* cov, int* supp_qual, int* supp, int* ref, int* max_qual)
{
	

	char reg_str[1000]; 
	sprintf(reg_str, "%s:%i-%i", v.chr, v.pos, v.pos+1); 

	//printf("samtools view %s  \n", fn_bam); 
	
	vector<const bam1_t*> reads;
	get_reads_from_bam(fn_bam, reg_str, &reads); 
	(*supp) = 0; 
	(*ref) = 0; 
	(*supp_qual) = 0; 
	(*max_qual) = 0; 
	for (uint32_t i=0; i<reads.size(); i++)
	{
		//std::cerr << "Current position " << reg_str << "\n";
		//uint8_t* rqual = bam1_qual(reads[i]); 
		//printf("%u  ", reads[i]->core.qual); 
		bool supp_ref; 
		bool supports = read_supports_variant(reads[i], v, &supp_ref);
		if (supports)
		{
			(*supp)++; 
			if (*(supp_qual)<reads[i]->core.qual)
				(*supp_qual) = reads[i]->core.qual; 
		}
		if (*(max_qual)<reads[i]->core.qual)
			(*max_qual) = reads[i]->core.qual; 

		if (supp_ref)
			(*ref)++; 

		//printf("read ID: %s %i, start_pos: %i\n", bam1_qname(reads[i]), supports, reads[i]->core.pos+1); 
	}
	//std::cerr << "Does read support variant? " << *supp << "\n";
	*cov = reads.size(); 
}


void split_variants(vector<variant>* vars)
{
	// split variants, where the ref and alt have the same length > 1 into 
	// individual SNVs
	int maxpos = 0; 
	for (uint32_t i=0; i<vars->size(); i++)
		if (vars->at(i).pos>maxpos) maxpos = vars->at(i).pos; 

	int cnt_new = 0; 
	int chrY = 0; 

	for (uint32_t i=0; i<vars->size(); i++)
	{
		if (strcmp(vars->at(i).chr, "chrY")==0)
		{
			chrY++; 
			vars->at(i).pos = maxpos+1000; // filter out
		}

		int len1 = strlen(vars->at(i).ref); 
		int len2 = strlen(vars->at(i).alt); 
		if (len1>1 && len1==len2)
		{
			for (int j=0; j<len1; j++)
			{
				if (vars->at(i).ref[j]!=vars->at(i).alt[j])
				{
					variant v = copy_var(vars->at(i)); 
					v.pos = vars->at(i).pos+j; 
					sprintf(v.ref, "%c", vars->at(i).ref[j]); 
					sprintf(v.alt, "%c", vars->at(i).alt[j]); 
					vars->push_back(v); 
					cnt_new++; 
				}
			}
			vars->at(i).pos = maxpos+1000; // sort them to the end of the array, then delete
		}
	}

	printf("[%s] sort %lu variants\n", __func__, vars->size()); 
	sort(vars->begin(), vars->end(), variant_compare_seq); 

	for (uint32_t i=0; i<vars->size(); i++)
	{
		if (vars->at(i).pos==maxpos+1000)
		{
			printf("[%s] replace %lu variants with %i new SNVs (chrY: %i)\n", __func__, vars->size()-i-chrY, cnt_new, chrY); 
			vars->resize(i); 
			break; 
		}
	}
}

bool comp_alt_allele(const variant& v1, const variant& v2)
        {
                if (strcmp(v1.chr, v2.chr)!=0)
                        return lexicographical_compare(v1.chr, v1.chr+strlen(v1.chr), v2.chr, v2.chr+strlen(v2.chr));
                if (v1.pos<v2.pos)
                        return true;
                if (v1.pos>v2.pos)
                        return false;
		if (strcmp(v1.ref, v2.ref)!=0)
			return lexicographical_compare(v1.ref, v1.ref+strlen(v1.ref), v2.ref, v2.ref+strlen(v2.ref));
	
                return lexicographical_compare(v1.alt, v1.alt+strlen(v1.alt), v2.alt, v2.alt+strlen(v2.alt));
        }

void merge_identical_variants(vector<variant>* vars) // changed .freq to .qual because for the ground truth variants, Jonas used the .qual field for the frequency, since they do not have a quality score anyways
{
	// after split_variants, some SNVs are listed several times in the ground truth SNV list with different frequencies
	// here they will be written as one SNV with the frequencies added

	vector<variant> unique_vars/* = vector<variant>(vars.size())*/; 
	int cnt_double=0;
	int cnt_unique=0;

	// sort vars according to chr, pos, and then alternative allele
	sort(vars->begin(),vars->end(),	comp_alt_allele);

	///*for (vector<variant>::const_iterator it = vars->begin(); it != vars->end(); ++it)
	//{
	//	std::cerr << it->chr << '\t' << it->pos << '\t' << it->ref << '\t' << it->alt << '\n';
	//}*/

	// now go through all of them
	variant last_variant_in_list = copy_var(vars->at(0)); // the first variant will be written here
	float currFreq=last_variant_in_list.qual; // and its frequency
	for (uint32_t i=1; i<vars->size(); ++i){ // start the list at the second variant
		bool isEqual=false;
		if(last_variant_in_list.pos==vars->at(i).pos){ // pos identical
			if(strcmp(last_variant_in_list.chr, vars->at(i).chr)==0){ // chr identical
				if(strcmp(last_variant_in_list.alt, vars->at(i).alt)==0){ // alternative allele identical
					//std::cerr << "[" << __func__ << "]"  << "These two variants are identical:\n";
					//std::cerr << "[" << __func__ << "]"  << last_variant_in_list.chr << '\t' << last_variant_in_list.pos << '\t' << last_variant_in_list.ref << '\t' << last_variant_in_list.alt << '\n';
					//std::cerr << "[" << __func__ << "]"  << vars->at(i).chr << '\t' << vars->at(i).pos << '\t' << vars->at(i).ref << '\t' << vars->at(i).alt << '\n';
					isEqual=true;
				}
			}
		}
	
		if(!isEqual){
			last_variant_in_list.qual=currFreq;
			unique_vars.push_back(last_variant_in_list);
			cnt_unique++;
			currFreq=vars->at(i).qual;
		} else {
			cnt_double++;
			currFreq+=vars->at(i).qual;
		}
		last_variant_in_list = copy_var(vars->at(i));
	}
	// after this loop, I still have to write the last variant into the vector
	last_variant_in_list.qual=currFreq;
	unique_vars.push_back(last_variant_in_list);
	// check which counter needs to be increased for this last one:
	uint32_t lastElement = vars->size()-1;
	if(currFreq==vars->at(lastElement).qual){ // the i is the last element of the vector vars
		cnt_unique++; // then we know nothing got added up
	} else {
		cnt_double++; // then we know the currFreq is more, because it is the sum of this and some identical variants before
	}

	printf("[%s] There are %d SNVs occuring double (positionwise)\n",__func__,cnt_double);
	printf("[%s] And there are %d SNVs occuring uniquely (positionwise)\n",__func__, cnt_unique);
	
	unique_vars.resize(cnt_unique);
	
	printf("[%s] Now we have in total %lu unique variants\n", __func__ , unique_vars.size()); 
	swap(*vars, unique_vars);
	for (unsigned int i=0; i<vars->size(); i++)
	{
		char* info = new char[1000]; // this will leak, fix if this is a limitation
		sprintf(info, "SRF=1,SRR=1,SAF=1,SAR=1,FREQ=%.10f", vars->at(i).qual); 
		vars->at(i).info = info; 
	}
	//sort(vars->begin(),vars->end(),	comp_alt_allele);
	sort(vars->begin(), vars->end(), variant_compare_seq); // sort it again Jonas' way such that it is conform with the remaining part of the script
}

void split_alt_alleles(vector<variant>* vars) // this is for the tools' vcf in a case like this: chr6	1612265	.	T	C,G
{
	// in case there are several alternate alleles listed, (separated by comma)
	// we split them up into multiple individual SNVs, at the same locus
	
	vector<variant> noCommaVars;
	noCommaVars.reserve(2 * vars->size()); // reserve the memory

	int cnt_vars_onlyComma=0, cnt_vars_toBeSplitted=0;
	for (uint32_t i=0; i<vars->size(); i++)
	{
		char *wholeAlt = vars->at(i).alt;
		if(strcmp(wholeAlt,",")==0){
			cnt_vars_onlyComma++;
			continue; // if the alternate allele is a comma, skip this variant; this should not occur at all
		}
		if( strchr(wholeAlt,',')!=NULL ) { // the alt contains a comma, which means that there are several alternate alleles for which we will create separate variants
			cnt_vars_toBeSplitted++;
			char *pointerToAlts;
			pointerToAlts = strtok (vars->at(i).alt,","); // get a pointer to the first element, when splitting it up according to commas
			while (pointerToAlts != NULL){ // this loops over all elements that were separated by comma
				char *currAlt = pointerToAlts;
				//printf ("%s\n",pointerToAlts);
				variant v = copy_var(vars->at(i)); 
				sprintf(v.alt, "%s", currAlt);
				noCommaVars.push_back(v);

				pointerToAlts = strtok (NULL, ",");
			}
		} else {
			noCommaVars.push_back(vars->at(i));
		}
	}
	printf("[%s] Out of %lu total variants, we found %d with only comma, and %d that contained a comma and were split.\n", __func__, vars->size(), cnt_vars_onlyComma, cnt_vars_toBeSplitted); 
	std::cout << "This resulted in a total of " << noCommaVars.size() << " variants\n";

	swap(*vars, noCommaVars);
	sort(vars->begin(), vars->end(), variant_compare_seq); // sort it again Jonas' way such that it is conform with the remaining part of the script
}


//double area_under_curve(double* xy, int len, bool reversed)
//{   
//
//	double area = 0.0;
//	
//	if (len==0)
//		return 0.0; 
//	
//	assert(xy); 
//
//	if (!reversed)
//	{   
//		for (int i=1; i<len; i++)
//			area += 0.5*(xy[2*i]-xy[2*(i-1)])*(xy[2*i+1]+xy[2*(i-1)+1]);
//	}
//	else
//	{   
//		for (int i=1; i<len; i++)
//			area += 0.5*(xy[2*i+1]-xy[2*(i-1)+1])*(xy[2*i]+xy[2*(i-1)]);
//	}
//
//	return area;
//}

void choose_high_pr(vector<pair<float, int> >* PR, vector<pair<float, int> >* result, float pr_cutoff)
{
	// sample from the prediction to compute a robust cutoff
	vector<int> cutoffs; 

	for (int iter=0; iter<50; iter++) 
	{
		vector<pair<float, int> > sample_PR; 
		for (size_t j=0; j<PR->size(); j++)
		{
			size_t idx = rand()%PR->size(); 
			sample_PR.push_back(PR->at(idx));
		}
		sort(sample_PR.begin(), sample_PR.end(), pair_comp);  

		int zeros=0; 
		int ones=0; 
		for (size_t i=0; i<sample_PR.size(); i++)
		{
			zeros += 1-sample_PR[i].second; 
			ones += sample_PR[i].second; 
		}

		size_t i; 
		for (i=0; i<sample_PR.size(); i++)
		{
			float frac = ((float) ones)/(zeros+ones); 
			if (frac>pr_cutoff)
			{
				break; 
			}
			zeros -= 1-sample_PR[i].second; 
			ones -= sample_PR[i].second; 
		}

		cutoffs.push_back(i); 

		//printf("cutoff: %lu (%lu) var: %f\n", i, PR->size(), emp_var(&cutoffs)); 
	}

	int cut_idx = my_prctile(&cutoffs, 0.5, false); // sort 

	printf("remove the lower %i (std: %.2f) variants (out of %lu -> %.2f%%)\n", cut_idx, sqrt(emp_var(&cutoffs)), PR->size(), ((float) cut_idx)/PR->size()*100); 
	
	for (size_t i=cut_idx; i<PR->size(); i++)
		result->push_back(PR->at(i)); 
}


float compute_auPRC_WriteCurve(vector<pair<float, int> >* PR, int total_ones, char* fn_PRcurve, bool writeCurve)
{
	//if (writeCurve){
	FILE* fd_PRcurve = fopen(fn_PRcurve, "w");
	assert(fd_PRcurve);
	//}

	printf("Total number of real positives that we count (without complex labels), i.e. true pos + false neg = %u\n",total_ones);
	// total_ones is the number of SNVs that are truly SNVs (excluding those labelled as complex) and that ideally
	// should be detected by the callers

	float tp = 0;
	double* curve = new double[2*PR->size()]; 
	//for(unsigned int i=0; i<PR->size(); i++)
	unsigned int cnt = 0;
	for(int i=PR->size()-1; i>=0; i--) // Jochen and Ariane noticed this loop started previously with the worst ranked SNVs; now we reverse the order
	{
		tp += PR->at(i).second==1;
		// precision 
		//curve[2*i] = tp/(i+1); 
		curve[2*cnt] = tp/(cnt+1); // Since the loop goes reverse, this should not be i anymore then; but cnt, which counts the number of SNVs looked at
		// recall
		//curve[2*i+1] = tp/total_ones; 
		curve[2*cnt+1] = tp/total_ones; // Since the loop goes reverse, this should not be i anymore then; but cnt, which counts the number of SNVs looked at
		if(writeCurve){
			//fprintf(fd_PRcurve, "%f\t%f\n", curve[2*i], curve[2*i+1]);
			fprintf(fd_PRcurve, "%f\t%f\n", curve[2*cnt], curve[2*cnt+1]); // Since the loop goes reverse, this should not be i anymore then; but cnt, which counts the number of SNVs looked at
		}
		if(i==0){
			break;
		}
		cnt++;
	}
	//if(writeCurve){
		fclose(fd_PRcurve);
	//}
	bool reversed = true; 
	double auPRC;
	if(PR->size()==0){
		auPRC=0.0;
	} else {
		auPRC = area_under_curve(curve, PR->size(), reversed); 
	}
	delete[] curve; 
	return auPRC; 
}



//float median(vector<float>* vec)
//{
//	sort(vec->begin(), vec->end()); 
//	float p = ((float) vec->size()-1)/2; 
//	if (vec->size()>ceil(p))
//		return (vec->at(floor(p))+vec->at(ceil(p)))/2; 
//	else
//		assert(false); 
//}

float get_precision_for_freq(vector<float>* uniq_freq, vector<pair<float, int> >* PR_freq, float freq, int* num)
{
	float lb=0.0, ub=0.0;
	for (unsigned int j=0; j<uniq_freq->size(); j++)
	{
		if (uniq_freq->at(j)==freq && j>0 && j<uniq_freq->size()-1)
		{
			lb = exp((log(uniq_freq->at(j))+log(uniq_freq->at(j-1)))/2); 
			ub = exp((log(uniq_freq->at(j))+log(uniq_freq->at(j+1)))/2); 
		}
		else if (uniq_freq->at(j)==freq && j==0)
		{
			lb = 0.0; 
			ub = exp((log(uniq_freq->at(j))+log(uniq_freq->at(j+1)))/2); 
		}
		else if (uniq_freq->at(j)==freq && j==uniq_freq->size()-1)
		{
			lb = exp((log(uniq_freq->at(j))+log(uniq_freq->at(j-1)))/2); 
			ub = 1.0; 
		}
	}
	//printf("freq: %f, lb:%f, ub:%f\n", freq, lb, ub); 
	*num=0; 
	int zeros_freq = 5; 
	int ones_freq = 0; 
	for (unsigned int j=0; j<PR_freq->size(); j++)
	{
		if (PR_freq->at(j).first>lb && PR_freq->at(j).first<=ub)
		{
			zeros_freq += 1-PR_freq->at(j).second; 
			ones_freq += PR_freq->at(j).second; 
			(*num)++; 
		}
	}
	float fPR = 0.0; 
	if ((zeros_freq+ones_freq)>0)
		fPR = ((float) ones_freq)/(zeros_freq+ones_freq); 
	return fPR; 
}


//void freq_distrib(vector<variant>* var)
//{
//	vector<float> freq;
//	for (unsigned int j=0; j<var->size(); j++)
//	{
//		freq.push_back(var->at(j).freq); 
//	}
//	assert(freq.size()>0); 
//	sort(freq.begin(), freq.end()); 
//	for (float p = 0; p<1; p+=0.05)
//	{
//		printf("%.2f:%.3f ", p, freq[round(p*(freq.size()-1))]); 
//	}
//	printf("\n"); 
//}

vector<float> read_abundance(char* fn_abundance)
{
	FILE* fd = fopen(fn_abundance, "r"); 
	assert(fd); 
	char line[10000]; 
	fgets(line, 10000, fd); 
	fclose(fd); 
	vector<char*> fields = my_str_tok(line, " ");

	float sum = 0.0; 
	vector<float> ret; 
	for (unsigned int i=0; i<fields.size(); i++)
	{
		ret.push_back(atof(fields[i])); 
		sum += ret.back(); 
	}
	for (unsigned int i=0; i<fields.size(); i++)
	{
		printf("%.8f, ", ret[i]/sum); 
		ret[i]/=(sum*2); // devide by two because we have this number for each of the two chromosomes 
	}
	printf("\n"); 
	return ret; 
}

//void combine_tree_vcf_files(char* vcf_dir, char* vcf_out, char* fn_abundance, char* normal_tag, vector<variant>* all_var)
void combine_tree_vcf_files(char* vcf_dir, char* fn_abundance, char* normal_tag, vector<variant>* all_var)
{
	vector<float> freq = read_abundance(fn_abundance); 
	assert(freq.size()==8); 

	for (int i=0; i<8; i++)
	{
		for (int j=0; j<=1; j++)
		{
			char fn_vcf[1000]; 
			sprintf(fn_vcf, "%s/%i_%i_TU_final.vcf", vcf_dir, i+7, j); 

			printf("read variants from file %s\n", fn_vcf);
			vector<variant> var;
			read_bcf(fn_vcf, &var, 0.0, 0.0, 0, normal_tag); 

			printf("sort new variants\n"); 
			sort(var.begin(), var.end(), variant_compare_seq); 	

			unsigned int p1 = 0; 
			unsigned int p2 = 0; 
			unsigned int num_all_var = all_var->size(); 
			while (p1<var.size())
			{
				if (p2==num_all_var)
				{
					var[p1].freq = freq[i];
					all_var->push_back(var[p1]); 
					p1++; 
				}
				else if (variant_equal_seq(var[p1], all_var->at(p2)))
				{
					assert(freq[i]+all_var->at(p2).freq<=1.01); 
					all_var->at(p2).freq += freq[i];
					if (all_var->at(p2).freq>1.0) all_var->at(p2).freq = 1.0; 
					p1++; 
				}
				else if (variant_compare_seq(var[p1], all_var->at(p2)))
				{
					// var[p1] is smaller than all_var->at(p2)
					var[p1].freq = freq[i];
					all_var->push_back(var[p1]); 
					p1++; 
				}
				else if (variant_compare_seq(all_var->at(p2), var[p1]))
				{
					p2++; 
				}
				else
				{
					assert(false); 
				}
			}
			printf("sort %lu variants\n", all_var->size()); 
			sort(all_var->begin(), all_var->end(), variant_compare_seq); 	
		}
	}
	for (unsigned int i=0; i<all_var->size(); i++)
	{
		char* info = new char[1000]; // this will leak, fix if this is a limitation
		sprintf(info, "SRF=1,SRR=1,SAF=1,SAR=1,FREQ=%.10f", all_var->at(i).freq); 
		all_var->at(i).qual = all_var->at(i).freq; 
		all_var->at(i).info = info; 
	}
	//write_vcf(vcf_out, all_var); // will be only written after splitting multinucleotide variants and then merging identical ones
}


void filter_variants(map<string, int*>* coverage_map, vector<variant>* var)
{
	int res = 10; 
	
	vector<variant> tmp; 
	map<string, int*>::iterator it; 
	for (unsigned int i=0; i<var->size(); i++)
	{
		it = coverage_map->find(string(var->at(i).chr)); 
		if (it == coverage_map->end())
		{
			printf("could not find field %s in map\n", var->at(i).chr);  
			exit(-1); 
		}
		if (it->second[var->at(i).pos/res]>=1)
			tmp.push_back(var->at(i)); 
	}
	printf("%lu out of %lu (%.2f%%) variants are in bed regions\n", tmp.size(), var->size(), ((float) 100*tmp.size())/(var->size())); 

	(*var) = tmp; 
}

void filter_variants_bam(map<string, vector<uint32_t> >* coverage_map, vector<variant>* var, int thresh)
{
	int res = 1; 
	
	vector<variant> tmp; 
	map<string, vector<uint32_t> >::iterator it; 
	for (unsigned int i=0; i<var->size(); i++)
	{
		it = coverage_map->find(string(var->at(i).chr)); 
		if (it == coverage_map->end())
		{
			printf("could not find field %s in map\n", var->at(i).chr);  
			exit(-1); 
		}
		if (((int) it->second[var->at(i).pos/res])>=thresh)
			tmp.push_back(var->at(i)); 
	}
	printf("%lu out of %lu (%.2f%%) variants are in regions with coverage >= %i\n", tmp.size(), var->size(), ((float) 100*tmp.size())/(var->size()), thresh); 

	(*var) = tmp; 
}


void parse_bed(map<string, int*>* coverage_map, const char* bamPath, const char* bedPath)
{
	char bed_ending[1000];
	sprintf(bed_ending, "/S04380219/S04380219_Regions.bed");
	//const char[] bed_ending = "/S04380219/S04380219_Regions.bed";
	char* fn_bed = (char *) std::malloc(1 + strlen(bedPath)+ strlen(bed_ending));
	strcpy(fn_bed, bedPath);
	strcat(fn_bed, bed_ending);
	std::cout  << "fn_bed = " << fn_bed << "\n";

	char bam_ending[1000];
	sprintf(bam_ending,"/alignments_sn_k1/bam/NO_final.10perc.bam");
	char* fn_bam = (char *) std::malloc(1 + strlen(bamPath)+ strlen(bam_ending));
	strcpy(fn_bam, bamPath);
	strcat(fn_bam,bam_ending); 
	std::cout << "fn_bam = " << fn_bam << "\n";

	bamFile fd1 = bam_open(fn_bam, "r");
	if (fd1==0)
	{   
		fprintf(stderr, "[%s] Could not open bam file: %s", __func__, fn_bam);
		exit(-1);
	}
	bam_header_t* header = bam_header_read(fd1);


	FILE* fd = fopen(fn_bed, "r"); 
	assert(fd); 
	char line[1000000]; 

	// get rid of header lines
	fgets(line, 1000000, fd); 
	assert(fgets(line, 1000000, fd)); 

	// story every k-th position only
	int res = 10; 

	for (int i=0; i<header->n_targets; i++)
	{   
		int len = header->target_len[i];
		char* chr = header->target_name[i];

		int* map = new int[len/res];
		memset(map, 0, len/res*sizeof(int));
		coverage_map->insert(pair<string, int*>(string(chr), map)); 
	}

	int cnt = 0; 
	while (fgets(line, 1000000, fd))
	{
		vector<char*> fields = my_str_tok(line, "\t"); 
		assert(fields.size()>=3); 
		char* chr = fields[0]; 
		int pos1 = atoi(fields[1]);  
		int pos2 = atoi(fields[2]);  

		map<string, int*>::iterator it; 
		it = coverage_map->find(string(chr)); 
		if (it == coverage_map->end())
		{
			printf("could not find field %s in map\n", chr);  
			exit(-1); 
		}
		cnt++; 
		int* map = it->second; 
		for (int i=pos1/res; i<pos2/res; i++)
		{
			map[i]++; 
		}
	} 
	printf("parsed %i regions ... ", cnt); 
	bam_header_destroy(header);
	bam_close(fd1);
	fclose(fd); 

}
bool var_compatible(variant v1, variant v2)
{
	if (strcmp(v1.chr, v2.chr)!=0)
		return false; 

	if (((int) (v1.pos+strlen(v1.ref)))<v2.pos && ((int) (v2.pos+strlen(v2.ref)))<v1.pos)
		return false; 

	if (v1.pos>v2.pos)
	{
		return var_compatible(v2, v1); 
	}


	vector<char*> vec_v1alt = my_str_tok(v1.alt, ","); 
	vector<char*> vec_v2alt = my_str_tok(v2.alt, ","); 

	bool match = false; 
	for (unsigned int i=0; i<vec_v1alt.size(); i++)
	{
		for (unsigned int j=0; j<vec_v2alt.size(); j++)
		{
			match = true; 

			char* v1alt = vec_v1alt[i]; 
			char* v2alt = vec_v2alt[j]; 
			int len1 = strlen(v1alt);
			int len2 = strlen(v2alt);
			
			if (v1.pos==v2.pos)
			{
				if (strcmp(v1alt, v2alt)==0)
					return true; 
				for (int i=0; i<len1 && i<len2; i++)
				{
					//if (v1alt[i] != v2alt[i])
					if (!IUB_compare(v1alt[i], v2alt[i]))
						match = false; 
				}
			}
			else if (v1.pos<v2.pos)
			{
				if (v2.pos+(len2-1)-v1.pos>=len1)
				{
					match = false; 
					continue; 
				}

				for (int i=0; i<len2; i++)
				{
					int l_pos = v2.pos-v1.pos+i; 
					assert(l_pos<len1); 
					//if (v1alt[l_pos]!=v2alt[i])
					if (!IUB_compare(v1alt[i], v2alt[i]))
					{
						match = false; 
					}
				}
			}
			if (match)
				return true; 
		}
	}
	return false;
}

bool is_complex_var(variant v)
{
	return strlen(v.ref)!=strlen(v.alt); 
}

bool check_nearby_indel(variant v, vector<variant>* all_var, int p2)
{
	unsigned int idx = p2+1; 
	int num_vars = 0; 
	while (idx<all_var->size())
	{	
		variant v2 = all_var->at(idx); 
		if (is_complex_var(v2))
			return true; 

		if (v2.pos>v.pos + 10)
			break; 

		num_vars += strlen(all_var->at(idx).alt); 
		idx++; 
	}
	idx = p2; 
	if (idx>0 && idx==all_var->size())
		idx--; 

	while (idx>0)
	{	
		variant v2 = all_var->at(idx); 

		if (is_complex_var(v2))
			return true; 

		if (v.pos>v2.pos + 10)
			break; 
		
		num_vars += strlen(all_var->at(idx).alt); 
		idx--; 
	}
	if (num_vars>=5)
		return true; 
	
	return false; 
}


bool full_var_compare(variant v, vector<variant>* all_var, int p2, int* found_idx)
{
	unsigned int idx = p2+1; 
	while (idx<all_var->size())
	{	
		variant v2 = all_var->at(idx); 

		if (var_compatible(v2, v))
		{
			*found_idx = idx; 
			return true; 
		}

		if (v2.pos>v.pos + 10)
			break; 
		idx++; 
	}
	idx = p2; 
	while (idx>0)
	{	
		variant v2 = all_var->at(idx); 

		if (var_compatible(v2, v))
		{
			*found_idx = idx; 
			return true; 
		}
		
		if (v.pos>v2.pos + 10)
			break; 
		idx--; 
	}
	
	return false; 
}

//evaluate_variants(fn_SN, fn_PR, &var, &var_tumor, &detected, fn_bam, fn_bam_normal,fn_SN_fdr01);
//float evaluate_variants(char* fn_SN, char* fn_PR, vector<variant>* var, vector<variant>* var_label, vector<int>* detected, char* fn_bam, char* fn_bam_normal,char* fn_SN_fdr01)

std::tuple<float, float, float>  evaluate_variants(char* fn_SN, char* fn_PR, vector<variant>* var, vector<variant>* var_label, vector<int>* detected, char* fn_bam, char* fn_bam_normal, char* fn_SN_fdr01, bool also_fp_fn)
{ // var_label are all ground truth tumor vars


	time_t timer;

	printf("sort %lu variants\n", var->size()); 
	sort(var->begin(), var->end(), variant_compare_seq); 	

	unsigned int p1 = 0; 
	unsigned int p2 = 0; 
	unsigned int num_var_label = var_label->size(); 
	
	unsigned int true_positive = 0; 
	unsigned int false_positive = 0; 
	unsigned int false_negative = 0; 

	vector<unsigned int> found(var_label->size(), 0);  // for each ground truth variable, get the information, whether it was found (=1), or not (=0)
	vector<unsigned int> all_found_idx(var->size(), 0); 
	vector<bool> near_indel(var->size(), 0); 

	// filter masks
	vector<unsigned int> label_is_complex(var_label->size(), 0); 
	vector<unsigned int> label_is_low_cov(var_label->size(), 0);
	vector<int> pred_hits_complex_label(var->size(), 0); 

	int found_idx = -1; 
	vector<pair<float, int> > SN; 
	vector<pair<float, int> > SN_qual; 
	vector<pair<float, int> > PR;  // for each detected variant, write here the quality score and the integer 0 or 1 that indicates whether it is a true pos or false pos
	vector<pair<float, int> > PR_freq;  // for each detected variant, write here the frequency and the integer 0 or 1 that indicates whether it is a true pos or false pos

	// check which ground truth variant is a 'complex' variant, i.e. imbalanced substitution
	// also if mode is also_fp_fn, check which ground truth variant is in a low coverage region
	for (unsigned int p2=0; p2<num_var_label; p2++)
	{
		label_is_complex[p2] = is_complex_var(var_label->at(p2)); 
		if(also_fp_fn){
			variant v = var_label->at(p2);
			int cov, qual, supp, ref, max_qual, cov_n, qual_n, supp_n, ref_n, max_qual_n;
			check_bam_qual(fn_bam, v, &cov, &qual, &supp, &ref, &max_qual);
			check_bam_qual(fn_bam_normal, v, &cov_n, &qual_n, &supp_n, &ref_n, &max_qual_n);
			if(cov < 25 || cov_n < 25)
				label_is_low_cov[p2] = 1;
		}
	}



	// check if variants are equivalent not just by comparing the two alternative sequences
	// here, the variants from the caller are checked whether they are true or false positives
	while (p1<var->size() || p2<num_var_label) // var is the vector with all detected variants, num_var_label the number of ground truth variants
	{
		if (p1==var->size()) // if p1 is at the end of counting, only increase p2 further
		{
			p2++;
		}
		else if (p2==num_var_label) // if p2 is at the end of counting, only increase p1
		{
			near_indel[p1] = check_nearby_indel(var->at(p1), var_label, p2); // returns true if variant close to imbalanced substitution or more than 4 SNVs
			PR.push_back(pair<float, int>(var->at(p1).qual, 0));
			PR_freq.push_back(pair<float, int>(var->at(p1).freq, 0));
			false_positive++;
			p1++; 
		}
		else if (variant_equal_seq(var->at(p1), var_label->at(p2)))
		{
			//if (var->at(p1).pos==183703238) assert(false); 
			PR.push_back(pair<float, int>(var->at(p1).qual, 1));
			PR_freq.push_back(pair<float, int>(var->at(p1).freq, 1));
			found[p2] = 1; 
			all_found_idx[PR.size()-1] = p2; 
			if (label_is_complex[p2] && pred_hits_complex_label[p1]!=-1)
				pred_hits_complex_label[p1] = 1; 
			else if (!label_is_complex[p2])
				pred_hits_complex_label[p1] =-1;

			p2++;
			p1++; 
		}
		else if (variant_compare_seq(var_label->at(p2), var->at(p1)))
		{
			p2++; 
		}
		//else if (var_label->at(p2).pos==var->at(p1).pos && var_label->at(p2).ref[0]!=var->at(p1).ref[0])
		//{
		//	// wrong ref base
		//	// this is a common problem of deepSNV, I contacted Moritz, how to deal with this. 
		//	// skip for now. 

		//	PR.push_back(pair<float, int>(var->at(p1).qual, 1));
		//	PR_freq.push_back(pair<float, int>(var->at(p1).freq, 1));
		//	found[p2] = 1; 
		//	all_found_idx[PR.size()-1] = p2; 

		//	p1++;
		//}
		else if (full_var_compare(var->at(p1), var_label, p2, &found_idx))
		{
			//if (var->at(p1).pos-var_label->at(found_idx).pos==0)
			//{
			//	printf("same pos variant match somehow:\n"); 
			//	printf("prediction:\t"); 
			//	print_vcf(stdout, var->at(p1)); 
			//	printf("label:     \t"); 
			//	print_vcf(stdout, var_label->at(found_idx)); 
			//} 
			//for (int i=0; i<var->at(p1).pos-var_label->at(found_idx).pos; i++)
			//	printf(" "); 
			//printf("%i %s (ref:%s)\n", var->at(p1).pos, var->at(p1).alt, var->at(p1).ref); 
			//for (int i=0; i<var_label->at(found_idx).pos-var->at(p1).pos; i++)
			//	printf(" "); 
			//printf("%i %s (ref:%s)\n", var_label->at(found_idx).pos, var_label->at(found_idx).alt,  var_label->at(found_idx).ref); 
			//printf("\n"); 
			//if (var->at(p1).pos==65564074) assert(false); 
			PR.push_back(pair<float, int>(var->at(p1).qual, 1));
			PR_freq.push_back(pair<float, int>(var->at(p1).freq, 1));
			found[found_idx] = 1; 
			all_found_idx[PR.size()-1] = found_idx; 

			if (label_is_complex[found_idx] && pred_hits_complex_label[p1]!=-1)
				pred_hits_complex_label[p1] = 1; 
			else if (!label_is_complex[found_idx])
				pred_hits_complex_label[p1] =-1;

			p1++; 
		}
		else if (variant_compare_seq(var->at(p1), var_label->at(p2)))
		{

			near_indel[p1] = check_nearby_indel(var->at(p1), var_label, p2); 
			PR.push_back(pair<float, int>(var->at(p1).qual, 0));
			PR_freq.push_back(pair<float, int>(var->at(p1).freq, 0));
			false_positive++; 
			p1++; 
		}
		else
		{
			assert(false);  // to make sure that there are never more cases than these considered here
		}
	}


	// create false negatives pie chart
	if (also_fp_fn)
	{
		char fn_FN_pie[1000];
		substitute(fn_SN, fn_FN_pie, "eval_SN", "eval_FN_pie_fix"); 
		char fn_FN_other_vcf[1000];
		substitute(fn_SN, fn_FN_other_vcf, "eval_SN", "FN_other.vcf"); 
		char fn_FN_low_cov_vcf[1000];
		substitute(fn_SN, fn_FN_low_cov_vcf, "eval_SN", "FN_low_cov.vcf"); 
		char fn_FN_low_qual_vcf[1000];
		substitute(fn_SN, fn_FN_low_qual_vcf, "eval_SN", "FN_low_qual.vcf"); 
		char fn_FN_near_indel_vcf[1000];
		substitute(fn_SN, fn_FN_near_indel_vcf, "eval_SN", "FN_near_indel.vcf"); 
		char fn_FN_low_support_vcf[1000];
		substitute(fn_SN, fn_FN_low_support_vcf, "eval_SN", "FN_low_support.vcf"); 
		char fn_FN_filtered_normal_vcf[1000];
		substitute(fn_SN, fn_FN_filtered_normal_vcf, "eval_SN", "FN_filtered_normal.vcf"); 
		char fn_FN_low_qual_normal_vcf[1000];
		substitute(fn_SN, fn_FN_low_qual_normal_vcf, "eval_SN", "FN_low_qual_normal.vcf"); 
		char fn_FN_low_cov_normal_vcf[1000];
		substitute(fn_SN, fn_FN_low_cov_normal_vcf, "eval_SN", "FN_low_cov_normal.vcf"); 
		char fn_FN_flag[1000];
		substitute(fn_SN, fn_FN_flag, "eval_SN", "FN_flags.txt");
		char fn_FN_absolute_counts[1000];
                substitute(fn_SN, fn_FN_absolute_counts, "eval_SN", "eval_FN_absolute_counts");


		FILE* fd_FN_pie = fopen(fn_FN_pie, "w"); 
		assert(fd_FN_pie); 
		FILE* fd_FN_other_vcf = fopen(fn_FN_other_vcf, "w"); 
		assert(fd_FN_other_vcf); 
		FILE* fd_FN_low_cov_vcf = fopen(fn_FN_low_cov_vcf, "w"); 
		assert(fd_FN_low_cov_vcf); 
		FILE* fd_FN_low_qual_vcf = fopen(fn_FN_low_qual_vcf, "w"); 
		assert(fd_FN_low_qual_vcf); 
		FILE* fd_FN_near_indel_vcf = fopen(fn_FN_near_indel_vcf, "w"); 
		assert(fd_FN_near_indel_vcf); 
		FILE* fd_FN_low_support_vcf = fopen(fn_FN_low_support_vcf, "w"); 
		assert(fd_FN_low_support_vcf); 
		FILE* fd_FN_filtered_normal_vcf = fopen(fn_FN_filtered_normal_vcf, "w"); 
		assert(fd_FN_filtered_normal_vcf); 
		FILE* fd_FN_low_qual_normal_vcf = fopen(fn_FN_low_qual_normal_vcf, "w"); 
		assert(fd_FN_low_qual_normal_vcf); 
		FILE* fd_FN_low_cov_normal_vcf = fopen(fn_FN_low_cov_normal_vcf, "w"); 
		assert(fd_FN_low_cov_normal_vcf); 
		FILE* fd_FN_flag = fopen(fn_FN_flag,"w");
		assert(fd_FN_flag);
		FILE* fd_FN_absolute_counts = fopen(fn_FN_absolute_counts,"w");
		assert(fd_FN_absolute_counts);
	

		int low_freq = 0; 
		int low_cov = 1; 
		int low_qual = 2; 
		int near_indel = 3; 
		int low_support = 4; 
		int filtered_normal = 5; 
		int low_qual_normal = 6; 
		int low_cov_normal = 7; 
		int other = 8; 
		vector<int>	all_counts(other+1, 0); 
		vector<int>	total_counts_category(other+1,0); // get the total numbers for each category

		for (uint32_t i=0; i<var_label->size(); i++)
		{

			if (found[i])
			{
				continue; 
			}
			if (label_is_complex[i]){ // imbalanced substitutions should not be counted as false negatives since each predicted variant that overlaps one, is not counted anyways // i.e. they are ignored for the SN and the FP pie charts, so they should be ignored here as well
				continue;
			}
			variant v = var_label->at(i); 
			if (v.qual<=0.25)// this is the frequency
			{
				all_counts[low_freq]++; 
				total_counts_category[low_freq]++;
				continue; 
			}
			
			vector<int>	counts_tmp(other+1, 0); 
			vector<int> flags(other+1,0); // to check which categories occurr how often absolutely; so here, for each variant i have a binary vector indicating which error category is true

			int cov, qual, supp, ref, max_qual, cov_n, qual_n, supp_n, ref_n, max_qual_n; 
			check_bam_qual(fn_bam, v, &cov, &qual, &supp, &ref, &max_qual);  
			check_bam_qual(fn_bam_normal, v, &cov_n, &qual_n, &supp_n, &ref_n, &max_qual_n);  

			std::cerr << "Current FN variant:\n";
			std::cerr << v.chr << "\t" << v.pos << "\t.\t" << v.ref << "\t" << v.alt << "\n";
			
			std::cerr << "cov\tqual\tsupp\tref\tmax_qual\n";
			std::cerr << cov <<"\t"<< qual <<"\t"<< supp <<"\t"<< ref <<"\t"<< max_qual <<"\n";

			std::cerr << "cov_n\tqual_n\tsupp_n\tref_n\tmax_qual_n\n";
			std::cerr << cov_n <<"\t"<< qual_n <<"\t"<< supp_n <<"\t"<< ref_n <<"\t"<< max_qual_n <<"\n";


			if (check_nearby_indel(v, var_label, i)){
				flags[near_indel]=1;
				counts_tmp[near_indel]++; 
				total_counts_category[near_indel]++;
				print_vcf(fd_FN_near_indel_vcf, v);
			}
			if (supp>=3 && qual<31){
				flags[low_qual]=1;
				counts_tmp[low_qual]++; 
				total_counts_category[low_qual]++;
				print_vcf(fd_FN_low_qual_vcf, v);
			}
			// if (binom_cdf(4, cov, v.qual)>0.5){// I expect to see less than 5 reads; v.qual is the frequency!!!
			if(cov<25){ 
				flags[low_cov]=1;
				counts_tmp[low_cov]++; 
				total_counts_category[low_cov]++;
				print_vcf(fd_FN_low_cov_vcf, v);
			}
			//if (binom_cdf(4, cov, v.qual)<=0.05 && supp<4){ // v.qual is the frequency!!! // if there are less than four reads supporting the variant, but only seeing 4 or less reads when having coverage "cov" and a frequency of "v.qual", is very unlikely
			if (cov>=25 && binom_cdf(4, cov, v.qual)<=0.05 && supp<4){ // v.qual is the frequency!!! // if there are less than four reads supporting the variant, but only seeing 4 or less reads when having coverage "cov" and a frequency of "v.qual", is very unlikely
				flags[low_support]=1;
				counts_tmp[low_support]++; 
				total_counts_category[low_support]++;
				//print_vcf(fd_FN_low_support_vcf, v); 
				fprintf(fd_FN_low_support_vcf, "%s\t%i\t%i\t%s\t%s\t%.3f\t%i\t%s\t.\n", v.chr, v.pos, supp ,v.ref, v.alt, v.qual, cov ,v.info); 
			}
			if (supp_n>1){ // n indicates normal
				flags[filtered_normal]=1;
				counts_tmp[filtered_normal]++; 
				total_counts_category[filtered_normal]++;
				print_vcf(fd_FN_filtered_normal_vcf,v);
			}
			// if (max_qual_n<31){
			if (cov_n>=25 && max_qual_n<31){ // n indicates normal  // otherwise it gets also incremented for variants with little or no coverage, but in the paper we say "... alignments with low mapping quality ... although the total coverage would be high enough"
				flags[low_qual_normal]=1;
				counts_tmp[low_qual_normal]++; 
				total_counts_category[low_qual_normal]++;
				print_vcf(fd_FN_low_qual_normal_vcf, v);
			}
			if (cov_n<25){
				flags[low_cov_normal]=1;
				counts_tmp[low_cov_normal]++; 
				total_counts_category[low_cov_normal]++;
				print_vcf(fd_FN_low_cov_normal_vcf, v);
			}
			
			int sum = my_sum(&counts_tmp);  

			if (sum==0)
			{
				counts_tmp[other]++;
				flags[other]=1;
				all_counts[other]++;
				total_counts_category[other]++;
				print_vcf(fd_FN_other_vcf, v); 
			}
			else
			{
				// randomly choose one case
				int idx = rand()%sum; 
				int c=0;
				int assigned = 0; 
				for (int i=0; i<other; i++)
				{
					if (counts_tmp[i]>0)
					{
						if (idx==c)
						{
							all_counts[i]++; 
							assigned++; 
						}
						c++; 
					}
				}
				assert(assigned==1); 
			}
			fprintf(fd_FN_flag, "%s\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", v.chr, v.pos, v.ref, v.alt , flags[low_freq], flags[low_cov], flags[low_qual], flags[near_indel], flags[low_support], flags[filtered_normal], flags[low_qual_normal], flags[low_cov_normal], flags[other]); // here we write out each variant together with the flag information, e.g. which error categories are true
	
		}

		printf("False neg: \n"); 
		printf("freq:%i ", all_counts[low_freq]); 
		printf("cov:%i ", all_counts[low_cov]); 
		printf("qual:%i ", all_counts[low_qual]); 
		printf("indel:%i ", all_counts[near_indel]); 
		printf("supp:%i ", all_counts[low_support]); 
		printf("filt_norm:%i ", all_counts[filtered_normal]); 
		printf("qual_norm:%i ", all_counts[low_qual_normal]); 
		printf("cov_norm:%i ", all_counts[low_cov_normal]); 
		printf("other:%i ", all_counts[other]); 
		

		int total_num_FN=0;
		for (int i=0; i<other; i++){
			fprintf(fd_FN_pie, "%i\t",all_counts[i]);
			total_num_FN+=all_counts[i];
		}
		fprintf(fd_FN_pie, "%i\n", all_counts[other]);
		total_num_FN+=all_counts[other];
		printf(" [TOTAL_FN] Total number of false negatives:%i \n", total_num_FN);
		
		// also print out the total counts
		for(int i=0; i<other; i++)
			fprintf(fd_FN_absolute_counts, "%i\t", total_counts_category[i]);
		fprintf(fd_FN_absolute_counts, "%i\n", total_counts_category[other]);

		fclose(fd_FN_pie); 
		fclose(fd_FN_other_vcf);
		fclose(fd_FN_low_cov_vcf);
		fclose(fd_FN_low_qual_vcf);
		fclose(fd_FN_near_indel_vcf);
		fclose(fd_FN_low_support_vcf);
		fclose(fd_FN_filtered_normal_vcf);
		fclose(fd_FN_low_qual_normal_vcf);
		fclose(fd_FN_low_cov_normal_vcf);
		fclose(fd_FN_flag);
		fclose(fd_FN_absolute_counts);

	}

	if (also_fp_fn)
	{
		std::cout<<"FP pie\n";
		vector<float> fdr; 
		//fdr.push_back(0.80);
		fdr.push_back(0.95); 
		//fdr.push_back(0.98); 
		fdr.push_back(0.99); 
		//fdr.push_back(0.995); 


		// write pie chart for false positives
		char fn_FP_pie[1000];
		sprintf(fn_FP_pie, "%s", fn_SN); 
		char* p_SN = strstr(fn_FP_pie, "eval_SN"); 
		assert(p_SN); 
		sprintf(p_SN, "eval_FP_pie_fix"); 
		FILE* fd_FP_pie = fopen(fn_FP_pie, "w"); 
		assert(fd_FP_pie); 

		char fn_FP_var_near_indel_vcf[1000];
		substitute(fn_SN, fn_FP_var_near_indel_vcf, "eval_SN", "FP_var_near_indel.vcf");
		char fn_FP_in_normal_vcf[1000];
		substitute(fn_SN, fn_FP_in_normal_vcf, "eval_SN", "FP_in_normal.vcf");
		char fn_FP_not_in_true_vcf[1000];
		substitute(fn_SN, fn_FP_not_in_true_vcf, "eval_SN", "FP_not_in_true.vcf");
		char fn_FP_very_high_cov_vcf[1000];
		substitute(fn_SN, fn_FP_very_high_cov_vcf, "eval_SN", "FP_very_high_cov.vcf");
		char fn_FP_other_vcf[1000];
		substitute(fn_SN, fn_FP_other_vcf, "eval_SN", "FP_other.vcf");
		char fn_FP_low_qual_vcf[1000];
		substitute(fn_SN, fn_FP_low_qual_vcf, "eval_SN", "FP_low_qual.vcf");
		char fn_FP_low_qual_normal_vcf[1000];
		substitute(fn_SN, fn_FP_low_qual_normal_vcf, "eval_SN", "FP_low_qual_normal.vcf");
		char fn_FP_low_cov_normal_vcf[1000];
		substitute(fn_SN, fn_FP_low_cov_normal_vcf, "eval_SN", "FP_low_cov_normal.vcf");
		char fn_FP_low_cov_vcf[1000];
		substitute(fn_SN, fn_FP_low_cov_vcf, "eval_SN", "FP_low_cov.vcf");
		char fn_FP_flag[1000];
		substitute(fn_SN, fn_FP_flag, "eval_SN", "FP_flags.txt");
		char fn_FP_absolute_counts[1000];
                substitute(fn_SN, fn_FP_absolute_counts, "eval_SN", "eval_FP_absolute_counts"); 



		FILE* fd_FP_var_near_indel_vcf = fopen(fn_FP_var_near_indel_vcf, "w"); 
		assert(fd_FP_var_near_indel_vcf); 
		FILE* fd_FP_in_normal_vcf = fopen(fn_FP_in_normal_vcf, "w"); 
		assert(fd_FP_in_normal_vcf); 
		FILE* fd_FP_not_in_true_vcf = fopen(fn_FP_not_in_true_vcf, "w"); 
		assert(fd_FP_not_in_true_vcf); 
		FILE* fd_FP_very_high_cov_vcf = fopen(fn_FP_very_high_cov_vcf, "w"); 
		assert(fd_FP_very_high_cov_vcf); 
		FILE* fd_FP_other_vcf = fopen(fn_FP_other_vcf, "w"); 
		assert(fd_FP_other_vcf); 
		FILE* fd_FP_low_qual_vcf = fopen(fn_FP_low_qual_vcf, "w"); 
		assert(fd_FP_low_qual_vcf); 
		FILE* fd_FP_low_qual_normal_vcf = fopen(fn_FP_low_qual_normal_vcf, "w"); 
		assert(fd_FP_low_qual_normal_vcf); 
		FILE* fd_FP_low_cov_normal_vcf = fopen(fn_FP_low_cov_normal_vcf, "w"); 
		assert(fd_FP_low_cov_normal_vcf); 
		FILE* fd_FP_low_cov_vcf = fopen(fn_FP_low_cov_vcf, "w"); 
		assert(fd_FP_low_cov_vcf); 
		FILE* fd_FP_flag = fopen(fn_FP_flag,"w");
		assert(fd_FP_flag);
		FILE* fd_FP_absolute_counts = fopen(fn_FP_absolute_counts,"w");
		assert(fd_FP_absolute_counts);
		
		for (unsigned int f=0; f<fdr.size(); f++)
		{
			std::cout<<"Current fdr: ";
			std::cout<<fdr[f]<< "\n";

			vector<unsigned int> found_cp = found; 
			vector<pair<float, int> > PR_sort; 
			for (unsigned int i=0; i<PR.size(); i++)
				PR_sort.push_back(pair<float,int>(PR[i].first, i)); 
			sort(PR_sort.begin(), PR_sort.end(), pair_comp); 
			int zeros=0; 
			int ones=0;
			int discard=0; 
			// find the largest set of variants from the top of the list, 
			// such that the precision is still above a given threshold
			for (unsigned int i=0; i<PR_sort.size(); i++)
			{
				if (pred_hits_complex_label[PR_sort[i].second] == 1)
					continue; 
				zeros += 1-PR[PR_sort[i].second].second; 
				ones += PR[PR_sort[i].second].second; 
			}
			unsigned int i; 
			for (i=0; i<PR_sort.size(); i++)
			{

				if (pred_hits_complex_label[PR_sort[i].second] == 1)
					continue; 

				float frac = ((float) ones)/(zeros+ones); 
				if (frac>fdr[f])
				{
					break; 
				}
				discard+=found[all_found_idx[PR_sort[i].second]]; 
				//if (found[all_found_idx[PR_sort[i].second]]==0)
				//	fp--; 
				found_cp[all_found_idx[PR_sort[i].second]]=0; 
				zeros -= 1-PR[PR_sort[i].second].second;
				ones -= PR[PR_sort[i].second].second; 
			}

			std::cout<<"Found the largest set of variants from the top of the list, such that the precision is still above a given threshold\n";
			std::cout<<"Now check the sources of false positives for the FP pie:\n";

			int var_near_indel = 0; 
			int in_normal = 1; 
			int not_in_true = 2; 
			int very_high_cov = 3; 
			int low_qual = 4;
			int low_qual_normal = 5;
			int low_cov = 6;
			int low_cov_normal = 7;
			int other = 8;
 
			vector<int> all_counts(other+1, 0);
			vector<int> total_counts_category(other+1,0); // to get the total numbers for each category

			for (unsigned int j=i; j<PR_sort.size(); j++) // we start counting at i, because this is the cutoff that was decided just before in the loop
			{
				vector<int> flags(other+1,0); // to check which categories occurr how often absolutely; here, for each variant, a binary vector is introduced indicating which error category is true

				if (PR[PR_sort[j].second].second!=0) // if it is a true positive, we continue to the next one
				{
					continue; 
				}
				if (j%1000==0)
					printf("\rFP: %u (%lu)", j, PR_sort.size()); 

				int cov, qual, supp, ref, max_qual,cov_n, qual_n, supp_n, ref_n, max_qual_n; 
				variant v = var->at(PR_sort[j].second); 
				check_bam_qual(fn_bam, v, &cov, &qual, &supp, &ref, &max_qual);  
				check_bam_qual(fn_bam_normal, v, &cov_n, &qual_n, &supp_n, &ref_n, &max_qual_n);
			
				std::cerr << "Current FP variant:\n";
				std::cerr << v.chr << "\t" << v.pos << "\t.\t" << v.ref << "\t" << v.alt << "\n";
		
				std::cerr << "cov\tqual\tsupp\tref\tmax_qual\n";
				std::cerr << cov <<"\t"<< qual <<"\t"<< supp <<"\t"<< ref <<"\t"<< max_qual <<"\n";

				std::cerr << "cov_n\tqual_n\tsupp_n\tref_n\tmax_qual_n\n";
				std::cerr << cov_n <<"\t"<< qual_n <<"\t"<< supp_n <<"\t"<< ref_n <<"\t"<< max_qual_n <<"\n";

				vector<int> counts_tmp(other+1, 0); 
				if (supp>=3 && qual<31){
					counts_tmp[low_qual]++;
					flags[low_qual]=1;
					total_counts_category[low_qual]++;
					if(fdr[f]<0.995 && fdr[f]>0.98)
						print_vcf(fd_FP_low_qual_vcf, v);
				}
				//if (max_qual_n<31){ // n indicates normal
				if (cov_n>=25 && max_qual_n<31){ // otherwise it gets also incremented for variants with little or no coverage, but in the paper we say "... alignments with low mapping quality ... although the total coverage would be high enough"
					counts_tmp[low_qual_normal]++; 
					flags[low_qual_normal]=1;
					total_counts_category[low_qual_normal]++;
					if(fdr[f]<0.995 && fdr[f]>0.98)
						print_vcf(fd_FP_low_qual_normal_vcf, v);
				}
				if (cov_n<25){
					counts_tmp[low_cov_normal]++;
					flags[low_cov_normal]=1;
					total_counts_category[low_cov_normal]++;
					if(fdr[f]<0.995 && fdr[f]>0.98)
						print_vcf(fd_FP_low_cov_normal_vcf, v);
				}
				if (cov<25){
					flags[low_cov]=1;
					counts_tmp[low_cov]++;
					total_counts_category[low_cov]++;
					if(fdr[f]<0.995 && fdr[f]>0.98)
						print_vcf(fd_FP_low_cov_vcf, v);
				}
				if (near_indel[PR_sort[j].second]){
					flags[var_near_indel]=1;
					counts_tmp[var_near_indel]++; 
					total_counts_category[var_near_indel]++;
					if(fdr[f]<0.995 && fdr[f]>0.98)
						print_vcf(fd_FP_var_near_indel_vcf, v);
				}
				if (cov>200){
					counts_tmp[very_high_cov]++; 
					flags[very_high_cov]=1;
					total_counts_category[very_high_cov]++;
					if(fdr[f]<0.995 && fdr[f]>0.98)
						print_vcf(fd_FP_very_high_cov_vcf, v);
				}
				if ((v.source&NOT_IN_TRUE_ALIGN)>0){
					counts_tmp[not_in_true]++; 
					flags[not_in_true]=1;
					total_counts_category[not_in_true]++;
					if(fdr[f]<0.995 && fdr[f]>0.98)
						print_vcf(fd_FP_not_in_true_vcf, v);
				}
				if ((v.source&MISSED_IN_NORMAL)>0){
					counts_tmp[in_normal]++; 
					flags[in_normal]=1;
					total_counts_category[in_normal]++;
					if(fdr[f]<0.995 && fdr[f]>0.98)
						print_vcf(fd_FP_in_normal_vcf,v);
				}

				int sum = my_sum(&counts_tmp);  

				if (sum==0){
					counts_tmp[other]++; 
					all_counts[other]++;
					total_counts_category[other]++;
					flags[other]=1;
					if(fdr[f]<0.995 && fdr[f]>0.98)
						print_vcf(fd_FP_other_vcf,v);
				}
				else
				{
					// randomly choose one case
					int idx = rand()%sum; 
					int c=0;
					int assigned = 0; 
					for (int i=0; i<other; i++)
					{
						if (counts_tmp[i]>0)
						{
							if (idx==c)
							{
								all_counts[i]++; 
								assigned++; 
							}
							c++; 
						}
					}
					assert(assigned==1); 
				}

				fprintf(fd_FP_flag, "%s\t%i\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", v.chr, v.pos, v.ref, v.alt , flags[var_near_indel], flags[in_normal], flags[not_in_true], flags[very_high_cov], flags[other], flags[low_qual], flags[low_qual_normal], flags[low_cov], flags[low_cov_normal]); // here we write out each variant together with the flag information, e.g. which error categories are true


			}
			std::cout<<"Write the error sources into the FP pie out file:\n";
			//int fp = zeros; 
			//printf("fraction: %f, zeros: %i (variable reg: %i in normal: %i not in true: %i high_cov: %i het_non_ref: %i other%i), ones: %i discard: %i i: %i\n", fdr[f], zeros, all_counts[var_near_indel], all_counts[in_normal], all_counts[not_in_true], all_counts[very_high_cov], all_counts[het_non_ref], all_counts[other], ones, discard, i); 

			//fprintf(fd_FP_pie, "%f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", fdr[f], all_counts[var_near_indel], all_counts[in_normal], all_counts[not_in_true], all_counts[very_high_cov], all_counts[het_non_ref], all_counts[other],all_counts[low_qual],all_counts[low_qual_normal],all_counts[low_cov],all_counts[low_cov_normal]); 
			fprintf(fd_FP_pie, "%f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", fdr[f], all_counts[var_near_indel], all_counts[in_normal], all_counts[not_in_true], all_counts[very_high_cov], all_counts[other],all_counts[low_qual],all_counts[low_qual_normal],all_counts[low_cov],all_counts[low_cov_normal]); 

			// also write out the total counts
			fprintf(fd_FP_absolute_counts, "%f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n", fdr[f], total_counts_category[var_near_indel], total_counts_category[in_normal], total_counts_category[not_in_true], total_counts_category[very_high_cov], total_counts_category[other],total_counts_category[low_qual],total_counts_category[low_qual_normal],total_counts_category[low_cov],total_counts_category[low_cov_normal]);



			printf("high scoring false positives: %i %f\n", fdr[f]==0.995, fdr[f]); 
			//if (f==3)
			//{
			//	printf("high scoring false positives:\n"); 
			//	for (unsigned int j=i; j<PR_sort.size(); j++)
			//	{
			//		if (PR[PR_sort[j].second].second==0)
			//		{
			//			variant v = var->at(PR_sort[j].second); 
			//			//if (PR_sort[i].first>150)
			//			if (near_indel[PR_sort[j].second]==0 && PR_sort[j].first>230)
			//				print_vcf(stdout, v); 
			//		}
			//	}
			//}

		}
		fclose(fd_FP_pie); 
		fclose(fd_FP_var_near_indel_vcf);
		fclose(fd_FP_in_normal_vcf);
		fclose(fd_FP_not_in_true_vcf);
		fclose(fd_FP_very_high_cov_vcf);
		fclose(fd_FP_other_vcf);
		fclose(fd_FP_low_qual_vcf);
		fclose(fd_FP_low_qual_normal_vcf);
		fclose(fd_FP_low_cov_vcf);
		fclose(fd_FP_low_cov_normal_vcf);
		fclose(fd_FP_flag);
		fclose(fd_FP_absolute_counts);
	}

	// SN for several cutoffs:
	{
		std::cout<<"SN start\n";
		vector<float> fdr; 

		fdr.push_back(0.70);
		fdr.push_back(0.80);
		fdr.push_back(0.90);
		fdr.push_back(0.95);
		fdr.push_back(0.99);


		FILE* fd = fopen(fn_SN, "w"); 
		if (!fd)
		{
			printf("[%s] could not write file %s\n", __func__, fn_SN); 
			exit(1); 
		}


		
		char fn_SN_variants_below_fdr[1000];
		substitute(fn_SN, fn_SN_variants_below_fdr, "eval_SN", "eval_SN_fdr05.vcf");
		FILE* fd_SN_variants_below_fdr = fopen(fn_SN_variants_below_fdr, "w"); 
		assert(fd_SN_variants_below_fdr); 

		char fn_SN_variants_below_fdr01[1000];
		substitute(fn_SN, fn_SN_variants_below_fdr01, "eval_SN", "eval_SN_fdr01.vcf");
		FILE* fd_SN_variants_below_fdr01 = fopen(fn_SN_variants_below_fdr01, "w"); 
		assert(fd_SN_variants_below_fdr01); 

		char fn_SN_variants_below_fdr10[1000];
		substitute(fn_SN, fn_SN_variants_below_fdr10, "eval_SN", "eval_SN_fdr10.vcf");
		FILE* fd_SN_variants_below_fdr10 = fopen(fn_SN_variants_below_fdr10, "w"); 
		assert(fd_SN_variants_below_fdr10); 
		
		char fn_SN_summaryFreq[1000];
		substitute(fn_SN, fn_SN_summaryFreq, "eval_SN", "eval_SN_summaryFreq");
		std::ofstream fd_SN_summaryFreq(fn_SN_summaryFreq);

		for (unsigned int f=0; f<fdr.size(); f++)
		{
			std::cout<<"Current fdr: ";
			std::cout<<fdr[f]<< "\n";

			vector<unsigned int> found_cp = found; 
			vector<pair<float, int> > PR_sort; 
			for (unsigned int i=0; i<PR.size(); i++)
				PR_sort.push_back(pair<float,int>(PR[i].first, i)); 
			sort(PR_sort.begin(), PR_sort.end(), pair_comp); 
			int zeros=0; 
			int ones=0;
			int discard=0; 
			// find the largest set of variants from the top of the list, 
			// such that the precision is still above a given threshold

			// count in the variables zeros and ones, how many true and false positives are in the total list of predictions
			unsigned int cntPredThatOverlapComplexGroundTruth=0;
			unsigned int cntTotalPred=PR_sort.size();
			for (unsigned int i=0; i<PR_sort.size(); i++)
			{
				if (pred_hits_complex_label[PR_sort[i].second] == 1){
					cntPredThatOverlapComplexGroundTruth++;
					continue;
				}
				zeros += 1-PR[PR_sort[i].second].second; 
				ones += PR[PR_sort[i].second].second; 
			}
			float percentagePredComplex=((float)cntPredThatOverlapComplexGroundTruth/cntTotalPred)*100;
			printf("PredThatOverlapComplexGroundTruth:%i\tTotalNumPred:%i\tPercentage:%.2f\n", cntPredThatOverlapComplexGroundTruth, cntTotalPred, percentagePredComplex);	
			unsigned int i; 


			int i_threshold;
			for (i=0; i<PR_sort.size(); i++)
			{

				if (pred_hits_complex_label[PR_sort[i].second] == 1)
					continue; 

				float frac = ((float) ones)/(zeros+ones); 
				if (frac>fdr[f])
				{
					i_threshold=i;
					break;
				}			
				
				discard+=found[all_found_idx[PR_sort[i].second]]; 
				//if (found[all_found_idx[PR_sort[i].second]]==0)
				//	fp--; 
				found_cp[all_found_idx[PR_sort[i].second]]=0;  // found_cp is the indicator for found or not found, with the respective cutoff
				zeros -= 1-PR[PR_sort[i].second].second;
				ones -= PR[PR_sort[i].second].second; 
			}

			if ( fdr[f] > 0.90 && fdr[f] < 0.99 ){ // 0.9500  -  for each caller, get the list printed with the variants that are below a certain FDR
				for (i=i_threshold; i<PR_sort.size(); i++)
				{
					variant v = var->at(PR_sort[i].second);	
					print_vcf(fd_SN_variants_below_fdr, v);
				}
			} else if (fdr[f] > 0.95){ // 0.9900000
				for (i=i_threshold; i<PR_sort.size(); i++)
				{
					variant v = var->at(PR_sort[i].second);	
					print_vcf(fd_SN_variants_below_fdr01, v);
				}
			} else if ( fdr[f] > 0.80 && fdr[f] < 0.95 ) { // 0.90000
				for (i=i_threshold; i<PR_sort.size(); i++)
				{
					variant v = var->at(PR_sort[i].second);
					print_vcf(fd_SN_variants_below_fdr10, v);
				}
			}

			std::cout<<"Found the largest set of variants from the top of the list, such that the precision is still above a given threshold\n";
			int fp = zeros; 
			std::cout<<"Now the SN file is being computed:\n";

			true_positive = 0; 
			false_negative = 0; 

			SN.clear(); 
			for (unsigned int i=0; i<found_cp.size(); i++)
			{
				if (label_is_complex[i])
					continue; 

				if (found_cp[i]==1)
				{
					if (detected)
						detected->at(i)++; 
					true_positive++;
					SN.push_back(pair<float, int>(var_label->at(i).qual, 1)); // this is not the quality but the frequency
				}
				else
				{
					false_negative++; 
					SN.push_back(pair<float, int>(var_label->at(i).qual, 0)); 
				}
			}
			float total_SN_at_prec=((float)true_positive)/(true_positive+false_negative);
			printf("cutoff:%f TP:%i, FP:%i FN:%i\tSN:%f\n", fdr[f], true_positive, fp, false_negative, total_SN_at_prec); 

			sort(SN.begin(), SN.end(), pair_comp);

			// SN now contains the frequencies of the ground truth variants and the information whether it is found or not found

			{
				int zeros = 0; 
				int ones = 0; 
				float prev = -1.0; 
				int num = 0; // num counts how many ground truth variants are there with the specified frequency

				vector<float> uniq_freq; 
				for (unsigned int i=0; i<SN.size(); i++)
				{
					uniq_freq.push_back(SN[i].first); 
				}
				sort(uniq_freq.begin(), uniq_freq.end()); 
				vector<float>::iterator it; 
				it = unique(uniq_freq.begin(), uniq_freq.end()); 
				uniq_freq.resize(it-uniq_freq.begin());
				
				// uniq_freq is now the unique vector of all ground truth frequencies
				printf("We have %lu unique ground truth frequencies.\n", uniq_freq.size());
		
				//std::cerr << "Current fdr round: " << fdr[f] << "\n";
				//std::cerr << "We have " << uniq_freq.size() << " unique ground truth frequencies. These are:\n";
				//for (unsigned int i=0; i<uniq_freq.size(); ++i){
				//	std::cerr << uniq_freq[i] << '\n';
				//}

				int num2; 
				int cnt_all=0;
				int cnt_tp=0;
				for (unsigned int i=0; i<SN.size(); i++)
				{
					if (prev != SN[i].first)
					{
						if (i>0)
						{
							float fPR = get_precision_for_freq(&uniq_freq, &PR_freq, prev, &num2);
							fprintf(fd, "%f\t%f\t%f\t%i\t%i\n", prev, ((float) ones)/(zeros+ones), fPR, num, num2);
						}
						zeros = 0;
						ones = 0; 
						num = 0; 
					}
					assert(prev<=1.0); 
					zeros+=SN[i].second==0; 
					ones+=SN[i].second==1; 
					num++; 
					prev = SN[i].first; 
					cnt_all++;
					cnt_tp+=SN[i].second==1;
				}
				float fPR = get_precision_for_freq(&uniq_freq, &PR_freq, prev, &num2);
				fprintf(fd, "%f\t%f\t%f\t%i\t%i\n", prev, ((float) ones)/(zeros+ones), fPR, num, num2);
				printf("We went through %d ground truth variants, of which %d were true positives.\n", cnt_all, cnt_tp);
			}
			fprintf(fd, "\n");

			// now the same thing as just before (writing out the sensitivity for a given frequency), but this time, for a smoother profile, we 
			// take a range of frequencies together, i.e. 0.05 range
			{
				vector<float> freq_breaks;
				//for(int i=1;i<=101;i++){
				//	freq_breaks.push_back((float)i/100);
				//}
				for(int i=1;i<=21;i++){
					freq_breaks.push_back((float)i/20);
				}
				//for(int i=1;i<=11;i++){
				//	freq_breaks.push_back((float)i/10);
				//}

				//std::cerr << "These are the frequency breaks at which the sensitivity will be evaluated:\n";
				//for(unsigned int i=0; i<freq_breaks.size(); ++i){
				//	std::cerr << freq_breaks[i] << "\n";
				//	fd_SN_summaryFreq << freq_breaks[i] << "\n";
				//}
				
				int zeros = 0; 
				int ones = 0; 
				int num = 0; // num counts how many ground truth variants are there with the specified frequency

				int cnt_all=0;
				int cnt_tp=0;
				int idx_freq=0;
				for (unsigned int i=0; i<SN.size(); i++)
				{
					if(SN[i].first >= freq_breaks[idx_freq]){ 
						//fprintf(fd_SN_summaryFreq, "%f\t%f\t%i\n", freq_breaks[idx_freq], ((float) ones)/(zeros+ones), num);
						fd_SN_summaryFreq << freq_breaks[idx_freq] << "\t" << ((float) ones)/(zeros+ones)  << "\t" << num  << std::endl;
						idx_freq++;
						zeros = 0;
						ones = 0;
						num = 0;
					}
					zeros+=SN[i].second==0;
					ones+=SN[i].second==1;
					num++;
					cnt_tp+=SN[i].second==1;
					cnt_all++;
				}
				//fprintf(fd_SN_summaryFreq, "%f\t%f\t%i\n", freq_breaks[idx_freq], ((float) ones)/(zeros+ones), num);
				fd_SN_summaryFreq << freq_breaks[idx_freq] << "\t" << ((float) ones)/(zeros+ones)  << "\t" << num  << std::endl;
				printf("We went through %d ground truth variants, of which %d were true positives.\n", cnt_all, cnt_tp);
			}
			fd_SN_summaryFreq << std::endl;

		}
		fclose(fd); 
		fclose(fd_SN_variants_below_fdr);
		fclose(fd_SN_variants_below_fdr01);
		fclose(fd_SN_variants_below_fdr10);

		fd_SN_summaryFreq.close();
	}


	// SN for when allowing up to 0.00 precision (i.e. taking all variants)
	{
		vector<float> fdr; 
		fdr.push_back(0.00);
		FILE* fd = fopen(fn_SN_fdr01, "w"); 
		if (!fd)
		{
			printf("[%s] could not write file %s\n", __func__, fn_SN_fdr01); 
			exit(1); 
		}

		for (unsigned int f=0; f<fdr.size(); f++)
		{
			std::cout<<"Current fdr: ";
			std::cout<<fdr[f]<< "\n";
			time(&timer);
			std::cout<<"Current time: "<<timer<<"\n";

			vector<unsigned int> found_cp = found; 
			vector<pair<float, int> > PR_sort; 
			for (unsigned int i=0; i<PR.size(); i++)
				PR_sort.push_back(pair<float,int>(PR[i].first, i)); 
			sort(PR_sort.begin(), PR_sort.end(), pair_comp); 
			int zeros=0; 
			int ones=0;
			int discard=0; 
			// find the largest set of variants from the top of the list, 
			// such that the precision is still above a given threshold
			for (unsigned int i=0; i<PR_sort.size(); i++)
			{
				if (pred_hits_complex_label[PR_sort[i].second] == 1)
					continue; 
				zeros += 1-PR[PR_sort[i].second].second; 
				ones += PR[PR_sort[i].second].second; 
			}

			unsigned int i; 
			for (i=0; i<PR_sort.size(); i++)
			{

				if (pred_hits_complex_label[PR_sort[i].second] == 1)
					continue; 

				float frac = ((float) ones)/(zeros+ones); 
				if (frac>fdr[f])
				{
					break; 
				}
				discard+=found[all_found_idx[PR_sort[i].second]]; 
				//if (found[all_found_idx[PR_sort[i].second]]==0)
				//	fp--; 
				found_cp[all_found_idx[PR_sort[i].second]]=0; 
				zeros -= 1-PR[PR_sort[i].second].second;
				ones -= PR[PR_sort[i].second].second; 
			}

			std::cout<<"Found the largest set of variants from the top of the list, such that the precision is still above a given threshold\n";
			int fp = zeros;
			std::cout<<"Now the SN file is being computed:\n";

			true_positive = 0; 
			false_negative = 0; 

			SN.clear(); 
			for (unsigned int i=0; i<found_cp.size(); i++)
			{
				if (label_is_complex[i])
					continue; 

				if (found_cp[i]==1)
				{
					if (detected)
						detected->at(i)++; 
					true_positive++;
					SN.push_back(pair<float, int>(var_label->at(i).qual, 1)); // this is not the quality but the frequency
				}
				else
				{
					false_negative++; 
					SN.push_back(pair<float, int>(var_label->at(i).qual, 0)); 
				}
			}
			printf("cutoff:%f TP:%i, FP:%i FN:%i\n", fdr[f], true_positive, fp, false_negative); 
			float cutoff_precision = ((float)true_positive) / (true_positive + fp);
			float cutoff_recall = ((float)true_positive) / (true_positive + false_negative);
			float cutoff_f1_score = ((float)2*true_positive) / (2*true_positive + false_negative + fp);
			printf("cutoff:%f\tprecision:%f\trecall:%f\tf1_score:%.2f\n", fdr[f], cutoff_precision, cutoff_recall,cutoff_f1_score);

			sort(SN.begin(), SN.end(), pair_comp);

			{
				int zeros = 0; 
				int ones = 0; 
				float prev = -1.0; 
				int num = 0;

				vector<float> uniq_freq; 
				for (unsigned int i=0; i<SN.size(); i++)
				{
					uniq_freq.push_back(SN[i].first); 
				}
				sort(uniq_freq.begin(), uniq_freq.end()); 
				vector<float>::iterator it; 
				it = unique(uniq_freq.begin(), uniq_freq.end()); 
				uniq_freq.resize(it-uniq_freq.begin());

				int num2; 
				for (unsigned int i=0; i<SN.size(); i++)
				{
					if (prev != SN[i].first)
					{
						if (i>0)
						{
							float fPR = get_precision_for_freq(&uniq_freq, &PR_freq, prev, &num2);
							fprintf(fd, "%f\t%f\t%f\t%i\t%i\n", prev, ((float) ones)/(zeros+ones), fPR, num, num2);
						}
						zeros = 0;
						ones = 0; 
						num = 0; 
					}
					assert(prev<=1.0); 
					zeros+=SN[i].second==0; 
					ones+=SN[i].second==1; 
					num++; 
					prev = SN[i].first; 
				}
				float fPR = get_precision_for_freq(&uniq_freq, &PR_freq, prev, &num2);
				fprintf(fd, "%f\t%f\t%f\t%i\t%i\n", prev, ((float) ones)/(zeros+ones), fPR, num, num2);
			}
			fprintf(fd, "\n");
		}
		fclose(fd); 
	} // end of eval_SN_fdr01

	//////////////////////////////
	// SN when taking into account only ground truth variants where there is sufficient coverage in normal and tumor:
	if (also_fp_fn)
	{
		std::cout<<"SN start\n";
		vector<float> fdr; 
		fdr.push_back(0.70);
		fdr.push_back(0.80);
		fdr.push_back(0.90);
		fdr.push_back(0.95);
		fdr.push_back(0.99);

		char fn_SN_covered_regions[1000];
		substitute(fn_SN, fn_SN_covered_regions, "eval_SN", "eval_SN_coveredRegions");
		FILE* fd_SN_covered_regions = fopen(fn_SN_covered_regions, "w"); 
		assert(fd_SN_covered_regions);

		char fn_SN_summaryFreq_coveredRegions[1000];
		substitute(fn_SN, fn_SN_summaryFreq_coveredRegions, "eval_SN", "eval_SN_summaryFreq_coveredRegions");
		std::ofstream fd_SN_summaryFreq_coveredRegions(fn_SN_summaryFreq_coveredRegions);

		for (unsigned int f=0; f<fdr.size(); f++)
		{
			std::cout<<"Current fdr: ";
			std::cout<<fdr[f]<< "\n";

			vector<unsigned int> found_cp = found; 
			vector<pair<float, int> > PR_sort; 
			for (unsigned int i=0; i<PR.size(); i++)
				PR_sort.push_back(pair<float,int>(PR[i].first, i)); 
			sort(PR_sort.begin(), PR_sort.end(), pair_comp); 
			int zeros=0; 
			int ones=0;
			int discard=0; 
			// find the largest set of variants from the top of the list, 
			// such that the precision is still above a given threshold

			// count in the variables zeros and ones, how many true and false positives are in the total list of predictions
			for (unsigned int i=0; i<PR_sort.size(); i++)
			{
				if (pred_hits_complex_label[PR_sort[i].second] == 1){
					continue;
				}
				zeros += 1-PR[PR_sort[i].second].second; 
				ones += PR[PR_sort[i].second].second; 
			}
			unsigned int i; 

			int i_threshold;
			for (i=0; i<PR_sort.size(); i++)
			{

				if (pred_hits_complex_label[PR_sort[i].second] == 1)
					continue; 

				float frac = ((float) ones)/(zeros+ones); 
				if (frac>fdr[f])
				{
					i_threshold=i;
					break;
				}			
				
				discard+=found[all_found_idx[PR_sort[i].second]]; 
				//if (found[all_found_idx[PR_sort[i].second]]==0)
				//	fp--; 
				found_cp[all_found_idx[PR_sort[i].second]]=0;  // found_cp is the indicator for found or not found, with the respective cutoff
				zeros -= 1-PR[PR_sort[i].second].second;
				ones -= PR[PR_sort[i].second].second; 
			}

			std::cout<<"Found the largest set of variants from the top of the list, such that the precision is still above a given threshold\n";
			int fp = zeros; 
			std::cout<<"Now the SN file is being computed:\n";

			true_positive = 0; 
			false_negative = 0; 

			SN.clear(); 
			for (unsigned int i=0; i<found_cp.size(); i++)
			{
				if (label_is_complex[i])
					continue; 
				if(label_is_low_cov[i])
					continue;

				if (found_cp[i]==1)
				{
					if (detected)
						detected->at(i)++; 
					true_positive++;
					SN.push_back(pair<float, int>(var_label->at(i).qual, 1)); // this is not the quality but the frequency
				}
				else
				{
					false_negative++; 
					SN.push_back(pair<float, int>(var_label->at(i).qual, 0)); 
				}
			}
			float total_SN_at_prec=((float)true_positive)/(true_positive+false_negative);
			printf("cutoff:%f TP:%i, FP:%i FN:%i\tSN:%f\n", fdr[f], true_positive, fp, false_negative, total_SN_at_prec); 

			sort(SN.begin(), SN.end(), pair_comp);

			// SN now contains the frequencies of the ground truth variants and the information whether it is found or not found

			{
				int zeros = 0; 
				int ones = 0; 
				float prev = -1.0; 
				int num = 0; // num counts how many ground truth variants are there with the specified frequency

				vector<float> uniq_freq; 
				for (unsigned int i=0; i<SN.size(); i++)
				{
					uniq_freq.push_back(SN[i].first); 
				}
				sort(uniq_freq.begin(), uniq_freq.end()); 
				vector<float>::iterator it; 
				it = unique(uniq_freq.begin(), uniq_freq.end()); 
				uniq_freq.resize(it-uniq_freq.begin());
				
				// uniq_freq is now the unique vector of all ground truth frequencies
				printf("We have %lu unique ground truth frequencies.\n", uniq_freq.size());
		
				//std::cerr << "Current fdr round: " << fdr[f] << "\n";
				//std::cerr << "We have " << uniq_freq.size() << " unique ground truth frequencies. These are:\n";
				//for (unsigned int i=0; i<uniq_freq.size(); ++i){
				//	std::cerr << uniq_freq[i] << '\n';
				//}

				int num2; 
				int cnt_all=0;
				int cnt_tp=0;
				for (unsigned int i=0; i<SN.size(); i++)
				{
					if (prev != SN[i].first)
					{
						if (i>0)
						{
							float fPR = get_precision_for_freq(&uniq_freq, &PR_freq, prev, &num2);
							fprintf(fd_SN_covered_regions, "%f\t%f\t%f\t%i\t%i\n", prev, ((float) ones)/(zeros+ones), fPR, num, num2);
						}
						zeros = 0;
						ones = 0; 
						num = 0; 
					}
					assert(prev<=1.0); 
					zeros+=SN[i].second==0; 
					ones+=SN[i].second==1; 
					num++; 
					prev = SN[i].first; 
					cnt_all++;
					cnt_tp+=SN[i].second==1;
				}
				float fPR = get_precision_for_freq(&uniq_freq, &PR_freq, prev, &num2);
				fprintf(fd_SN_covered_regions, "%f\t%f\t%f\t%i\t%i\n", prev, ((float) ones)/(zeros+ones), fPR, num, num2);
				printf("We went through %d ground truth variants, of which %d were true positives.\n", cnt_all, cnt_tp);
			}
			fprintf(fd_SN_covered_regions, "\n");

			// now the same thing as just before (writing out the sensitivity for a given frequency), but this time, for a smoother profile, we 
			// take a range of frequencies together, i.e. 0.05 range
			{
				vector<float> freq_breaks;
				//for(int i=1;i<=101;i++){
				//	freq_breaks.push_back((float)i/100);
				//}
				for(int i=1;i<=21;i++){
					freq_breaks.push_back((float)i/20);
				}
				//for(int i=1;i<=11;i++){
				//	freq_breaks.push_back((float)i/10);
				//}

				//std::cerr << "These are the frequency breaks at which the sensitivity will be evaluated:\n";
				//for(unsigned int i=0; i<freq_breaks.size(); ++i){
				//	std::cerr << freq_breaks[i] << "\n";
				//	fd_SN_summaryFreq << freq_breaks[i] << "\n";
				//}
				
				int zeros = 0; 
				int ones = 0; 
				int num = 0; // num counts how many ground truth variants are there with the specified frequency

				int cnt_all=0;
				int cnt_tp=0;
				int idx_freq=0;
				for (unsigned int i=0; i<SN.size(); i++)
				{
					if(SN[i].first >= freq_breaks[idx_freq]){ 
						//fprintf(fd_SN_summaryFreq, "%f\t%f\t%i\n", freq_breaks[idx_freq], ((float) ones)/(zeros+ones), num);
						fd_SN_summaryFreq_coveredRegions << freq_breaks[idx_freq] << "\t" << ((float) ones)/(zeros+ones)  << "\t" << num  << std::endl;
						idx_freq++;
						zeros = 0;
						ones = 0;
						num = 0;
					}
					zeros+=SN[i].second==0;
					ones+=SN[i].second==1;
					num++;
					cnt_tp+=SN[i].second==1;
					cnt_all++;
				}
				//fprintf(fd_SN_summaryFreq, "%f\t%f\t%i\n", freq_breaks[idx_freq], ((float) ones)/(zeros+ones), num);
				fd_SN_summaryFreq_coveredRegions << freq_breaks[idx_freq] << "\t" << ((float) ones)/(zeros+ones)  << "\t" << num  << std::endl;
				printf("We went through %d ground truth variants, of which %d were true positives.\n", cnt_all, cnt_tp);
			}
			fd_SN_summaryFreq_coveredRegions << std::endl;

		}
		fclose(fd_SN_covered_regions); 
		fd_SN_summaryFreq_coveredRegions.close();
	}
	//////////////////////////////


	//float auPRC = 0.0;
	float auPRC05, auPRC10, auPRC100 = 0.0; 
	// auPRC05 for 5% FDR
	{
		sort(PR.begin(), PR.end(), pair_comp);
	
		int num_labels = 0; 
		for (size_t i=0; i<label_is_complex.size(); i++)
			num_labels+=1-label_is_complex[i]; 

		vector<pair<float, int> > top_k; 

		float myPrecisionCutoff = 0.95; // 5% FDR

		printf("The precision cutoff is %f.\n",myPrecisionCutoff);
		choose_high_pr(&PR, &top_k, myPrecisionCutoff); 

		char fn_PRcurve[1000];
		substitute(fn_SN, fn_PRcurve, "eval_SN", "PRcurve"); // to write out the precision recall curve for 100% FDR, i.e. whole list
		//auPRC = compute_auPRC_WriteCurve(&top_k, num_labels,fn_PRcurve); // this is to write out the whole precision recall curve

		//auPRC = compute_auPRC(&top_k, num_labels); 	// for the auPRC, only the top predictions are used to not penalize tools that only output little but with high precision
								// as opposed to tools that output a lot and therefore reach a higher sensitivity with lower precision
								// use this to compute the auPRC values for the plot with the coverage, 5%, 10%
		bool writeCurve = false; 
		auPRC05 = compute_auPRC_WriteCurve(&top_k, num_labels,fn_PRcurve,writeCurve);
		printf("auPRC: %f\n", auPRC05); 
		assert(!isnan(auPRC05)); 
		
		if (false) {
			FILE* fd_pr = fopen(fn_PR, "w"); 
			assert(fd_pr); 

			{
				vector<float> q_vals; 
				for (unsigned int i=0; i<PR.size(); i++)
				q_vals.push_back(PR[i].first); 

				vector<float> ptile = prctile<float>(&q_vals, 30);  

				int zeros = 0; 
				int ones = 0; 
				unsigned int prc_idx = 1; 

				for (unsigned int i=0; i<PR.size(); i++)
				{
					assert(prc_idx<ptile.size()); 
					if (PR[i].first > ptile[prc_idx])
					{
						if (i>0)
						{
							fprintf(fd_pr, "%f\t%f\n", ptile[prc_idx], ((float) ones)/(zeros+ones));
						}
						zeros = 0;
						ones = 0; 
						prc_idx++; 
					}
					zeros+=PR[i].second==0; 
					ones+=PR[i].second==1; 
				}
				fprintf(fd_pr, "%f\t%f\n", ptile.back(), ((float) ones)/(zeros+ones));
			}
			fclose(fd_pr); 
		}
	}

	// auPRC10 for 10% FDR
	{
		sort(PR.begin(), PR.end(), pair_comp);
	
		int num_labels = 0; 
		for (size_t i=0; i<label_is_complex.size(); i++)
			num_labels+=1-label_is_complex[i]; 

		vector<pair<float, int> > top_k; 

		float myPrecisionCutoff = 0.90; // 10% FDR

		printf("The precision cutoff is %f.\n",myPrecisionCutoff);
		choose_high_pr(&PR, &top_k, myPrecisionCutoff); 

		char fn_PRcurve[1000];
		substitute(fn_SN, fn_PRcurve, "eval_SN", "PRcurve"); // to write out the precision recall curve for 100% FDR, i.e. whole list
		//auPRC = compute_auPRC_WriteCurve(&top_k, num_labels,fn_PRcurve); // this is to write out the whole precision recall curve

		//auPRC = compute_auPRC(&top_k, num_labels); 	// for the auPRC, only the top predictions are used to not penalize tools that only output little but with high precision
								// as opposed to tools that output a lot and therefore reach a higher sensitivity with lower precision
								// use this to compute the auPRC values for the plot with the coverage, 5%, 10%
		bool writeCurve = false; 
		auPRC10 = compute_auPRC_WriteCurve(&top_k, num_labels,fn_PRcurve,writeCurve);
		printf("auPRC: %f\n", auPRC10); 
		assert(!isnan(auPRC10)); 	
	}
	
	// 100% FDR - i.e. the whole precision recall curve
	{
		sort(PR.begin(), PR.end(), pair_comp);
	
		int num_labels = 0; 
		for (size_t i=0; i<label_is_complex.size(); i++)
			num_labels+=1-label_is_complex[i]; 

		vector<pair<float, int> > top_k; 
		//choose_high_pr(&PR, &top_k, 0.99); 


		float myPrecisionCutoff = 0.00; //  all of it

		printf("The precision cutoff is %f.\n",myPrecisionCutoff);
		choose_high_pr(&PR, &top_k, myPrecisionCutoff); 

		char fn_PRcurve[1000];
		substitute(fn_SN, fn_PRcurve, "eval_SN", "PRcurve"); // to write out the precision recall curve
		bool writeCurve = true;
		auPRC100 = compute_auPRC_WriteCurve(&top_k, num_labels,fn_PRcurve,writeCurve); // this is to write out the whole precision recall curve

		printf("auPRC: %f\n", auPRC100); 
		assert(!isnan(auPRC100)); 
	}
	//return auPRC; 
	return std::make_tuple(auPRC05, auPRC10, auPRC100);
}


int is_indel(variant var)
{
	int len1 = strlen(var.ref); 
	int len2 = strlen(var.alt); 
	if (len1==1 && len2==1)
	{
		return 0; 
	}
	else
		return 1; 

	if (len1>len2)
		return len1-len2; 
	else
		return len2-len1; 
}


void filter_snv_indel(vector<variant>* var, int indel)
{
	vector<variant> tmp; 
	for (unsigned int i=0; i<var->size(); i++)
	{
		if (indel==0 && is_indel(var->at(i))==0)
			tmp.push_back(var->at(i)); 
		else if (indel>0)
		{
			int len_diff = is_indel(var->at(i)); 
			if (len_diff>0 && len_diff<=indel)
				tmp.push_back(var->at(i)); 
		}
	}
	printf("[%s] number of variants: %lu -> %lu\n", __func__, var->size(), tmp.size()); 
	*var = tmp; 
}

#endif
