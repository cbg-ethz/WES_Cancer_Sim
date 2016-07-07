# WES_Cancer_Sim

# sim_cancer:
- contains the scripts for the simulation, starting the variant callers, the combinations and overlaps, as well as the evaluation
- to run the analysis, refer to the 'README' in the respective subfolder in the following order:
	1. ./sim_cancer/simulate_WES_cancerSample/
	2. ./sim_cancer/var_callers/
	3. ./sim_cancer/combis_overlaps/
	4. ./sim_cancer/eval_vcf/

To do before running: add the paths to the paths.sh and paths.R file, and compile all programs

# cancer_data_pipeline:
- contains the scripts for alignment, variant calling, checking the coverage, and processing of the fastq/bam/variant files
- they are called by scripts in the sim_cancer subfolders


