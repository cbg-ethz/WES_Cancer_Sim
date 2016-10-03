# sim_cancer:

This folder contains the scripts for the simulation, starting the variant callers, the combinations and overlaps, as well as the evaluation.

To run the analysis, refer to the 'README' in the respective subfolder in the following order:

	1. ./simulate_WES_cancerSample/

	2. ./var_callers/

	3. ./combis_overlaps/

	4. ./eval_vcf/

To do before running: add the paths to the paths.sh and paths.R file, and compile all programs

# cancer_data_pipeline:
- contains the scripts for alignment, variant calling, checking the coverage, and processing of the fastq/bam/variant files
- they are called by scripts in the other subfolders


