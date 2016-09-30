###
# scripts to run for variant calling:
start_callers.sh: just prints out the commands to run for variant calling in default mode, variant calling with parameter tuning, repeated subsampling, local realignment, checking the coverage, ... ; prints the run commands for run_variant_callers.sh, tune_variant_callers.sh

###
# scripts to run for filtering and re-formatting of the output
filterOutGermlineComma_somSniper.sh: filters out germline mutations from the somaticSniper output; calls the program Extract_Somatic_Mult_Alleles.java

add_qual_in_VCF_for_ordering.sh: reformats the output of VarScan2 and somaticSniper. It writes the quality score of the variant into the 6th column, like is the case for a usual vcf file; calls the programs RewriteVarScanVCF.jar and RewriteSomaticSniperQualField.java

filterOutGermline_GATK_SAMtools.sh: filters out the germline mutations from the vcfs of GATK and SAMtools; calls the program WriteDifference_TU_NO_asSomaticFilter.java 

format_HPCaller.sh: re-formats the GATK HP vcf, because sometimes the column with the "GT:AD:DP:GQ:PL" details is malformed; calls the program FormatHPCallerGATK_TUonly.java

analysis_wCNVs_reformat.sh: all re-formats for the CNV analysis

###
# Now, the vcfs of the variant callers can be evaluated. see ../eval_vcf/

# or combined/intersected. see ../combis_overlaps/

