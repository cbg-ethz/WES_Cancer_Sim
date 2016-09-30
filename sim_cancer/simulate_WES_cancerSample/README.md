###
# order of scripts to run for simulation:
prepare_input_bams.sh: creates the two bam files from each of the tumor and normal bam files, one with all multimappers, and one without multimappers

create_cancer_vcfs_as_input.sh:	runs freebayes on all bam files

intersect_extract_merge_vcfs.sh: from the final tumor and normal vcf files, generate the intersections as desired; this script calls compute_overlap.cpp

run_create_cancer_vcfs.sh: inserts the variants into the tree and then, when you have a vcf for each clone, insert them into the reference genome, for each clone; this script calls create_cancer_vcfs.py and create_ref_from_vcf.cpp

run_wessim_all.sh: for each of the cancer clone genomes, it generates reads and aligns them; this script calls run_wessim.sh, process_alignments.sh, process_true_alignments.sh; first: make sure that 'select_true_alignment.cpp' in the folder ../eval_vcf/ is compiled


###
# Now, the variant callers can be started. see ../var_callers/

###
## these are additional scripts which can be run:

check_coverage.sh: for coverage statistics and qualimap

repeated_subsampling.sh: to assert the variability of the performance, repeated subsampling of the whole (100%) bam file



## Analysis with CNVs and aneuploidies

To add the call the CNVs, and add CNVs and aneuploidies to the simulation:

./analysis_wCNVs.sh

(uses the scripts generate_vcf_files_allSubclones.sh; run_control_freec_for_cnvs.sh; SNV_Overlapper_Ordered_PlusAlleles.java; distribute_cnvs_on_tree.py; extract_variants_in_cnv_regions_losses.sh; extract_variants_in_cnv_regions_gains.sh; combine_and_sort_subclones_variants_n_regions.sh; create_ref_from_vcf; run_wessim_wCNVs.sh; process_true_alignments_wCNVs.sh; process_alignments_wCNVs.sh)

Now, the variant callers can be started. see ../var_callers



