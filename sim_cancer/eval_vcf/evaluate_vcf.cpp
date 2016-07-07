#include <assert.h>
#include <tuple>
#include <string>
	using std::string;
#include <string.h>
#include <stdio.h>
#include <vector>
	using std::vector;
#include <map>
#include <unordered_set>
#include <unordered_map>
	using std::map;
	using std::unordered_multimap;
#include "gtf_tools.h"
#include <algorithm>
#include "region.h"
#include <limits> 
#include <zlib.h>
#include <math.h>

#include <ctype.h>

#include <bam.h>


#include "config.h"
#include "variant.h"
#include "bioinf_utils.h"
#include "bcf_utils.h"
#include "string_utils.h"
#include "file_utils.h"
#include "vector_math.h" 
#include "bam_tools.h"
#include "evaluate_vcf.h"
#include "stats_tests.h"
#include <iostream>
#include <fstream>


void usage()
{
	printf("Usage: evaluate_vcf --out-dir <dirname> --vcf-dir <dirname> --eval-dir <dirname> <fn_vcf_1> ...\n");  
}

int main(int argc, char* args[])
{
	//for (int i=0; i<=30; i++)
	//{
	//	printf("%i %f\n", i, binom_cdf(i, 30, 0.25)); 
	//}
	//exit (0); 

	char* normal_tag = NULL; 
	if (!normal_tag)
	{
		normal_tag = new char[1000]; 
		sprintf(normal_tag, "NO"); 
	}

	// params
	char out_dir[1000]; out_dir[0] = '\0';
	char vcf_dir[1000]; vcf_dir[0] = '\0';
	char eval_dir[1000]; eval_dir[0] = '\0';
	char bamPath[1000]; bamPath[0] = '\0';
	char bedPath[1000]; bedPath[0] = '\0';
	char* fn_bam = NULL; 
	char* fn_bam_normal = 0; 
	int indel = 0; // maximal size of insertions and deletions
	int max_cnt_normal = 0; 
	bool also_fp_fn=true;
	bool binom_test=false;

	vector<char*> vcf_files; 
	for (int i=1; i<argc; i++)
	{
		if (strcmp(args[i], "--out-dir")==0)
		{
			assert(argc>i+1); i++;
			sprintf(out_dir, "%s", args[i]); 
		}
		else if (strcmp(args[i], "--vcf-dir")==0)
		{
			assert(argc>i+1); i++;
			sprintf(vcf_dir, "%s", args[i]); 
		}
		else if (strcmp(args[i], "--eval-dir")==0)
		{
			assert(argc>i+1); i++;
			sprintf(eval_dir, "%s", args[i]); 
		}
		else if (strcmp(args[i], "--indel")==0)
		{
			assert(argc>i+1); i++;
			indel = atoi(args[i]); 
		}
		else if (strcmp(args[i], "--max-cnt-normal")==0)
		{
			assert(argc>i+1); i++;
			max_cnt_normal = atoi(args[i]);
		}
		else if (strcmp(args[i], "--bam")==0)
		{
			assert(argc>i+1); i++;
			fn_bam = args[i];
		}
		else if (strcmp(args[i], "--bam-normal")==0)
		{
			assert(argc>i+1); i++;
			fn_bam_normal = args[i];
		}
		else if(strcmp(args[i], "--only-sn-auprc")==0)
		{
			also_fp_fn=false; // since the FP and FN categories are time consuming, we only do this for the cont 20 and perc 50 cases
		}
		else if(strcmp(args[i], "--binom-test")==0)
		{
			binom_test=true; // filter the vcf list with predictions with a binomial test to remove potential germline mutations
		}
		else if (strcmp(args[i], "--pathToBed")==0)
		{
			assert(argc>i+1); i++;
			sprintf(bedPath, "%s", args[i]); 
		}
		else if(strcmp(args[i], "--pathToBam")==0)
		{
			assert(argc>i+1); i++;
			sprintf(bamPath, "%s", args[i]); 
		}
		else
			vcf_files.push_back(args[i]); 
	}
	if (strlen(out_dir)==0 || strlen(vcf_dir)==0 || strlen(eval_dir)==0)
	{
		usage(); 
		exit(0); 
	}

	//vector<variant> all_varxx; 
	//read_bcf(vcf_files[0], &all_varxx, 0.0, 0.0, 0, normal_tag);
	//exit(0); 

	printf("parse bed file ... "); 
	map<string, int*> coverage_map; 
	parse_bed(&coverage_map, bamPath, bedPath);
	printf("done\n"); 

	map<string, vector<uint32_t> > bam_cov_map; 
	char fn_abundance[1000]; 
	sprintf(fn_abundance, "%s/abundance.txt", out_dir); 
	char vcf_eval[1000]; 
	sprintf(vcf_eval, "%s/all_tumor_eval.vcf", out_dir);  // used to be written by the function "combine_tree_vcf_files" (without splitting or merging variants)
							      // now it is written only after splitting and merging
	char vcf_eval_bed[1000];
	char vcf_eval_normal[1000];
	sprintf(vcf_eval_bed, "%s/all_tumor_eval_bed_no_normal.vcf", out_dir); 
	sprintf(vcf_eval_normal, "%s/all_normal_eval.vcf", out_dir);	
	char vcf_eval_normal_bed[1000];
	sprintf(vcf_eval_normal_bed, "%s/all_normal_eval_bed.vcf", out_dir);

	char vcf_normal0[1000]; 
	sprintf(vcf_normal0, "%s/0_NO_final.vcf", vcf_dir); 
	char vcf_normal1[1000]; 
	sprintf(vcf_normal1, "%s/1_NO_final.vcf", vcf_dir); 
	

	vector<variant> var_tumor; 
	vector<variant> var_normal; 
	vector<variant> all_var; 


	printf("Load ground truth normal variants from %s and %s \n", vcf_normal0, vcf_normal1);
	// load normal ground truth
	read_bcf(vcf_normal0, &var_normal, 0.0, 0.0, 0, normal_tag);
	read_bcf(vcf_normal1, &var_normal, 0.0, 0.0, 0, normal_tag);
	//sort(var_normal.begin(), var_normal.end(), variant_compare_seq);
	// split variants and sort: 
	split_variants(&var_normal); 
	merge_identical_variants(&var_normal); 
	
	if(!fexist(vcf_eval_normal)){
		printf("Write the ground truth normal varinats into %s \n", vcf_eval_normal);
		write_vcf(vcf_eval_normal,&var_normal);	
	}

	if (fexist(vcf_eval_bed))
	{
		read_bcf(vcf_eval_bed, &var_tumor, 0.0, 0.0, 0, normal_tag);
		split_variants(&var_tumor); 
		merge_identical_variants(&var_tumor); 
	}
	else
	{
		printf("Create tumor vcf file! ");
		if (fexist(vcf_eval))
		{
			read_bcf(vcf_eval, &all_var, 0.0, 0.0, 0, normal_tag);
			split_variants(&all_var);
			merge_identical_variants(&all_var); 
		}
		else
		{
			//combine_tree_vcf_files(vcf_dir, vcf_eval, fn_abundance, normal_tag, &all_var); 
			combine_tree_vcf_files(vcf_dir, fn_abundance, normal_tag, &all_var);
			split_variants(&all_var); 
			merge_identical_variants(&all_var);
			write_vcf(vcf_eval, &all_var);
		}

		// load normal ground truth
		if (var_normal.size()==0)
		{
			printf("Load ground truth normal variants from %s and %s \n", vcf_normal0, vcf_normal1);
			read_bcf(vcf_normal0, &var_normal, 0.0, 0.0, 0, normal_tag);
			read_bcf(vcf_normal1, &var_normal, 0.0, 0.0, 0, normal_tag);
			//sort(var_normal.begin(), var_normal.end(), variant_compare_seq); 	
			split_variants(&var_normal); 
			merge_identical_variants(&var_normal); 
		}
		var_tumor = vector<variant>(all_var.size()+var_normal.size()); 

		vector<variant>::iterator it; 
		it = set_difference(all_var.begin(), all_var.end(), var_normal.begin(), var_normal.end(), var_tumor.begin(), variant_compare_seq); 
		var_tumor.resize(it-var_tumor.begin());
		printf("tumor%lu normal%lu set_diff%lu\n", all_var.size(), var_normal.size(), var_tumor.size()); 

		filter_variants(&coverage_map, &var_tumor); // coverage from bed regions

		write_vcf(vcf_eval_bed, &var_tumor); 
	}

	printf("filter snv (%i): label\n", indel);  

	vector<int> detected(var_tumor.size(), 0); 

	printf("evaluate %lu vcf files\n", vcf_files.size()); 
	fflush(stdout);
	for (unsigned int i=0; i<vcf_files.size(); i++)
	{
		vector<variant> var; 	
		printf("evaluate : %s\n", vcf_files[i]); 
		char fname[1000]; 
		char fn_SN[1000]; 
		char fn_SN_fdr01[1000];
		char fn_PR[1000]; 
		char fname4[1000]; 
		char fname4001[1000];
		char fname4000[1000];
		sprintf(fname, vcf_files[i]); 
		printf("basename: %s\n", basename(fname)); 
		sprintf(fn_SN, "%s/%s_indel%i.eval_SN", eval_dir, basename(fname), indel);      
		sprintf(fn_SN_fdr01, "%s/%s_indel%i.eval_SN_fdr01", eval_dir, basename(fname), indel);
		sprintf(fn_PR, "%s/%s_indel%i.eval_PR", eval_dir, basename(fname), indel);  
		sprintf(fname4001, "%s/%s_indel%i.eval_auPRC_10", eval_dir, basename(fname), indel);  
		sprintf(fname4, "%s/%s_indel%i.eval_auPRC_5", eval_dir, basename(fname), indel);
		sprintf(fname4000, "%s/%s_indel%i.eval_auPRC_100", eval_dir, basename(fname), indel);
		
		if (strstr(vcf_files[i], "samvar") && ! strstr(vcf_files[i], "samvar_1_2"))
		{
			char fn_normal[1000]; 
			sprintf(fn_normal, "%s", vcf_files[i]); 
			char* final = strstr(vcf_files[i], "final"); 
			char* TU = strstr(fn_normal, "TU."); 
			sprintf(TU, "NO_%s", final); 
			printf("fn_normal: %s\n", fn_normal); 
			assert(fexist(fn_normal)); 

			vector<variant> normal; 
			vector<variant> tumor; 
			read_bcf(fn_normal, &normal, 0.0, 0.0, 0, normal_tag);
			read_bcf(vcf_files[i], &tumor, 0.0, 0.0, 0, normal_tag);
			split_variants(&normal); 
			split_variants(&tumor); 
			//sort(normal.begin(), normal.end(), variant_compare_seq);
			//sort(tumor.begin(), tumor.end(), variant_compare_seq);

			var = vector<variant>(normal.size()+tumor.size()); 
			vector<variant>::iterator it; 
			it = set_difference(tumor.begin(), tumor.end(), normal.begin(), normal.end(), var.begin(), variant_compare_seq); 
			var.resize(it-var.begin());
		}
		else
		{
			read_bcf(vcf_files[i], &var, 0.0, 0.0, 0, normal_tag);
			split_alt_alleles(&var);
			//sort(var.begin(), var.end(), variant_compare_seq);
			split_variants(&var); 
			char vcf_predicted[1000];
			sprintf(vcf_predicted, "%s/%s_afterSplitting.vcf", eval_dir, basename(fname)); 
			write_vcf(vcf_predicted, &var);
		}

		// remove variants outside of bed regions
		filter_variants(&coverage_map, &var); 
		// remove variants of other types than specified by indel 
		printf("filter snv (%i): prediction\n", indel);  
		filter_snv_indel(&var, indel); 
		//filter_variants_bam(&bam_cov_map, &var_tumor, 50); 


		// if mode is binom_test, filter out mutations with the binomial filter that removes potential germline variants
		// if the coverage and variant count in the normal cannot be explained by a sequencing error, it is likely a germline mutation and is filtered out
		if(binom_test){
			printf("Now, the binomial test is performed and likely germline mutations are filtered out.\n");
			vector<variant> no_germline_vars; // no need to define size because with push_back, it is just reallocated on the fly
			no_germline_vars.reserve(var.size()); // here you tell it that it will have a size of at most this much, but reserving is not really necessary
			unsigned int count_somatic = 0;
			for(unsigned int var_idx=0; var_idx<var.size(); ++var_idx){
				//printf("This is variant number %i.\n",var_idx);
				const variant& v = var[var_idx]; // this, in contrast to "variant v = ...", does not copy the whole structure, but is just a nickname for the same thing, kind of like a pointer
				printf("Current variant: %s\t%i\n",v.chr, v.pos);
				int cov_n=0, qual_n=0, supp_n=0, ref_n=0, max_qual_n=0;
				check_bam_qual(fn_bam_normal, v, &cov_n, &qual_n, &supp_n, &ref_n, &max_qual_n);
				printf("cov_n = %i,\tsupp_n = %i\n",cov_n,supp_n);
				//printf("Bam file checked.\n");

				float var_rate = 0.005; // this is the sequencing error rate we assume
				/*estimated illumina error rate ~0.1% -> rate = 0.001
				here we use 0.005 to be more on the conservative side, i.e. the varCount needs to be higher in order to have a rejection = 
				germline classification = filter out variant
				NullHypothesis: varCount in normal is just sequencing error
				AltHypothesis: varCount in normal is higher and cannot be due to seq error -> must be germline
				in R: p.val=binom.test(varCount, coverage, p = 0.005, alternative = "greater",conf.level = 0.95)$p.val
				*/

				float pval = -1.0;
				if(cov_n == 0 || supp_n == 0){
					pval=1;
				} else if(cov_n == supp_n){
					pval=0;
				} else {
					pval = binom_cdf_alternative(supp_n, cov_n, var_rate,"greater"); // perform binomial test; null-hypothesis: only seq error; alt-hypothesis: germline mut
				}
				//printf("The p-value is %f.\n",pval);
				if(pval >= 0.05 ){ // Include as somatic 
					no_germline_vars.push_back(v);
					count_somatic++;
				}
				//if (pval < 0.05) { // if pval is less than 0.05, reject the null-hypoth. and classify as germline -> filter out later
			}
			//printf("These are the variants in no_germline_vars:\n");
			//for(unsigned int var_idx=0; var_idx<no_germline_vars.size();var_idx++){
			//	const variant& v = no_germline_vars[var_idx];
			//	printf("Current variant: %s\t%i\n",v.chr, v.pos);
			//}
			
			printf("Only keep somatic variants that passed: %i variants.\n",count_somatic);
			var.clear(); 
			vector<variant>::iterator it;
			unsigned int count_no_germline = 0;
			for (it=no_germline_vars.begin(); it!=no_germline_vars.end(); it++)
			{
					count_no_germline++;
					var.push_back(*it); 
					//const variant& v = *it;
					//printf("Current variant: %s\t%i\n",v.chr, v.pos);
			}
			assert(count_somatic==count_no_germline);
			char fn_vcf_pred_no_germline[1000];
			sprintf(fn_vcf_pred_no_germline, "%s/%s_binomTestSomatic.vcf", eval_dir, basename(fname));
			if(!fexist(fn_vcf_pred_no_germline)){
				printf("Write the binomial tested varinats into %s \n", fn_vcf_pred_no_germline);
				write_vcf(fn_vcf_pred_no_germline,&var);
			}
		}

		// mark variants `in normal` 
		if (also_fp_fn) // this is only necessary if the FN and FP categories are computed
		{
			vector<variant> var_no_normal = vector<variant>(var.size()+var_normal.size()); 

			vector<variant>::iterator it; 
			it = set_difference(var.begin(), var.end(), var_normal.begin(), var_normal.end(), var_no_normal.begin(), variant_compare_seq); 
			var_no_normal.resize(it-var_no_normal.begin());
			
			vector<variant> only_normal = vector<variant>(var_no_normal.size()+var.size()); 
			it = set_difference(var.begin(), var.end(), var_no_normal.begin(), var_no_normal.end(), only_normal.begin(), variant_compare_seq);
			only_normal.resize(it-only_normal.begin());

			printf("only_normal.size(): %lu, var_no_normal.size(): %lu, var.size(): %lu \n", only_normal.size(), var_no_normal.size(), var.size()); 

			// give those variants, that are found in the normal vcf a flag and reinsert then int the the prediction
			for (it=only_normal.begin(); it!=only_normal.end(); it++)
				it->source += MISSED_IN_NORMAL;

			var_no_normal.insert(var_no_normal.begin(), only_normal.begin(), only_normal.end()); 
			//printf("tumor%lu normal%lu set_diff%lu\n", var.size(), var_normal.size(), var_no_normal.size()); 

			sort(var_no_normal.begin(), var_no_normal.end(), variant_compare_seq);

			// filter out all variants, where the caller has observed confirmation in normal
			var.clear(); 
			for (it=var_no_normal.begin(); it!=var_no_normal.end(); it++)
			{
				//assert(it->conf_cnt_normal>=0); 
				//if (it->conf_cnt_normal <= max_cnt_normal) // we do not remove any mutations just because the vcf says that there are also reads in the normal sample that support this
					var.push_back(*it); 
			}
			printf(" remove %i (out of %lu) variants with nonzero evidence in normal\n", ((int) var_no_normal.size())-(int) var.size(), var_no_normal.size()); 
			//var = var_no_normal; 
		}



		if (also_fp_fn) // this is only necessary if the FN and FP categories are computed
		{ 
			if (strstr(vcf_files[i], "perc.") && ! strstr(vcf_files[i], "perc.true"))
			{
				// also load the variants detected with true alignments
				char fn_var_true[1000]; 
				sprintf(fn_var_true, "%s", vcf_files[i]); 
				char* p_perc = strstr(fn_var_true, "perc."); 
				assert(p_perc); 
				int num_written = sprintf(p_perc, "perc.true."); 
			
				char* p_perc_plus = strstr(vcf_files[i], "perc.") + 5; 
	
				sprintf(p_perc+num_written, "%s", p_perc_plus); 
				printf("load variants from fn_var_true: %s\n", fn_var_true); 
	
				vector<variant> var_true; 
				read_bcf(fn_var_true, &var_true, 0.0, 0.0, 0, normal_tag);
				split_alt_alleles(&var_true);
				//sort(var.begin(), var.end(), variant_compare_seq);
				split_variants(&var_true); 
				char vcf_predicted_true[1000];
				sprintf(vcf_predicted_true, "%s/%s_afterSplitting.true.vcf", eval_dir, basename(fname)); 
				write_vcf(vcf_predicted_true, &var_true);
	
	
				filter_variants(&coverage_map, &var_true); 
				// remove variants of other types than specified by indel 
				printf("true align: filter snv (%i): prediction\n", indel);  
				filter_snv_indel(&var_true, indel); 
	
				vector<variant> not_in_true(var.size()+var_true.size());
				vector<variant>::iterator it; 
				it = set_difference(var.begin(), var.end(), var_true.begin(), var_true.end(), not_in_true.begin(), variant_compare_seq); 
				not_in_true.resize(it-not_in_true.begin());
				printf("found %lu variants (out of %lu) (%.4f%%) which are not in the true variants set\n", not_in_true.size(), var.size(), ((float) not_in_true.size()/var.size())*100); 
	
				vector<variant> var_common(var.size()+var_true.size());
				it = set_intersection(var.begin(), var.end(), var_true.begin(), var_true.end(), var_common.begin(), variant_compare_seq); 
				var_common.resize(it-var_common.begin()); 
	
				// give variants a flag and reinsert into the list
				for (uint32_t j=0; j<not_in_true.size(); j++)
				{
					//if (not_in_true[j].qual > 200 && rand()>0.99 )
					//	print_vcf(stdout, not_in_true[j]); 
					not_in_true[j].source += NOT_IN_TRUE_ALIGN;  
				}
	
				var_common.insert(var_common.end(), not_in_true.begin(), not_in_true.end()); 
				sort(var_common.begin(), var_common.end(), variant_compare_seq);
	
				var = var_common; 
			}
		} 

		printf("Start the evaluate_variants function.\n");

		//float auPRC; 
		float auPRC05, auPRC10, auPRC100 = 0.0;
		std::tuple<float, float, float> auPRC;
		if (strstr(vcf_files[i], "50perc"))
			//auPRC = evaluate_variants(fn_SN, fn_PR, &var, &var_tumor, &detected, fn_bam, fn_bam_normal,fn_SN_fdr01); 
			auPRC = evaluate_variants(fn_SN, fn_PR, &var, &var_tumor, &detected, fn_bam, fn_bam_normal,fn_SN_fdr01,also_fp_fn);
		else
			//auPRC = evaluate_variants(fn_SN, fn_PR, &var, &var_tumor, NULL, fn_bam, fn_bam_normal, fn_SN_fdr01);
			auPRC = evaluate_variants(fn_SN, fn_PR, &var, &var_tumor, NULL, fn_bam, fn_bam_normal, fn_SN_fdr01,also_fp_fn);

		auPRC05 = std::get<0>(auPRC);
		auPRC10 = std::get<1>(auPRC);
		auPRC100 = std::get<2>(auPRC);

		{
		FILE* fd_prc = fopen(fname4, "w"); // 5% FDR
		assert(fd_prc); 
		fprintf(fd_prc, "%.10f\n", auPRC05); 
		fclose(fd_prc);
		}

		{
		// 10% FDR
		FILE* fd_prc = fopen(fname4001, "w"); // already done for all
		assert(fd_prc); 
		fprintf(fd_prc, "%.10f\n", auPRC10); 
		fclose(fd_prc); 
		}

		{
		// 100% FDR - meaning we take all variants
		FILE* fd_prc = fopen(fname4000, "w"); 
		assert(fd_prc); 
		fprintf(fd_prc, "%.10f\n", auPRC100); 
		fclose(fd_prc); 
		}
		
	}
	
	
	FILE* fd = fopen(paste("sss", eval_dir, "/", "missed.txt"), "w"); 
	assert(fd); 
	FILE* fd2 = fopen(paste("ss", eval_dir, "/found.txt"), "w"); 
	assert(fd2); 
	for (unsigned int i=0; i<var_tumor.size(); i++)
	{
		//fprintf(fd, "%i\t%.5f\t%s\t%i\n", detected[i], var_tumor[i].qual, var_tumor[i].chr, var_tumor[i].pos); 
		if (detected[i]==0)//&& var_tumor[i].qual==0.5)
			print_vcf(fd, var_tumor[i]); 
		else if (detected[i]>0)// && var_tumor[i].qual==0.5)
			print_vcf(fd2, var_tumor[i]); 
	}
	fclose(fd); 
	fclose(fd2); 

	exit(0); 
}












