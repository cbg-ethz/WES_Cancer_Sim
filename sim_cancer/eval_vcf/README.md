###
# order of scripts/commands to run for evaluation:
make: to compile evaluate_vcf.cpp select_true_alignment.cpp, and filter_somatic_binom_test_JointSNVMix2.cpp; first: change the path to samtools to your path to samtools

filter_somatic_binomTest_JSM2.sh: if you want to assess the germline filter, this script runs the germline filter 'binomial test' on the JointSNVMix2 vcf; the other variant callers get this germline filter in the main evaluation script evaluate_vcf.cpp; JointSNVMix2 gets this separate treatment, because it ouputs > 19 million variants, and it would take way too long to check; the bam file at all these loci --> we rely on the variant count numbers and coverages in the JointSNVMix2 vcf; this script starts filter_somatic_binom_test_JointSNVMix2.cpp

run_eval.sh: the main script to start all evaluations; first: set the runs you want to evaluate to 'true'; this script starts evaluate_vcf.cpp for the evaluation

###
# Now, the results can be plotted:
plot_SN_curves.R: plots the sensitivity as a function of variant allele frequency for a fixed precision

plot_FN_FP_absolute_counts.R: plots the error profiles of false negative and false positive error sources, also plots the correlation of error profiles

plot_vennDiagram.R:  plots the vennDiagrams, where the numbers have been generated with the scripts in ../combis_overlaps/

plot_auPRC_summaryBarplot.R: plots the summary bar plot for a direct comparison of combis, and individual tools

plot_auPRC.sh: plots the area under precision recall curve (auPRC) as a function of coverage and contamination; calls the script plot_auPRC.R

plot_PRC.sh: plots the precision recall curves for all tools; starts the scripts plot_PRC.R, and plot_PRC_altogether.R

plot_coverage_profiles.sh: plots the coverage profiles of the 20% contamination, 50% coverage file; starts plot_coverage_profiles.R

plot_auPRC_repeatedSubsampling.sh: plots the auPRC of all the repeated subsamplings in a boxplot; starts plot_auPRC_repeatedSubsampling.R

