#!/bin/bash

# Please note: The user has to have Seqan installed, to run the analysis
# In your shell: export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${gitDir}/sim_cancer/eval_vcf/imported_scripts/hdf5/hdf5-1.8.11/hdf5/lib/


# directories
gitDir=`pwd`
dir_=$gitDir
SAMDIR=<pathTo>/samtools-0.1.19
SNV_DIR=$gitDir
base_dir=<pathTo>/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60
vcf_sim_dir=<pathTo>/sim_cancer_exome/data/subsample/final_vcfs_0.25
pathToBed=<pathTo>/agilent/
pathToBam=$base_dir
overlapDir=<pathTo>
localToolsDir=<pathTo>


# tools
bowtie_index_tool=$localToolsDir/bowtie2-2.1.0/bowtie2-build
bowtie=$localToolsDir/bowtie2-2.1.0/bowtie2
picard_tools="java -Xmx2g -jar $localToolsDir/picard-tools-1.92/"
fastqc=$localToolsDir/fastqc-v0.10.1/FastQC/fastqc
coverageBed=$localToolsDir/BEDTools2.17.0/bedtools-2.17.0/bin/coverageBed
intersectBed=$localToolsDir/BEDTools2.17.0/bedtools-2.17.0/bin/intersectBed
somatic_sniper=$localToolsDir/somatic-sniper/somaticSniper_v1.0.4/somatic-sniper-1.0.4/build-common/bin/bam-somaticsniper
bcftools=$localToolsDir/samtools-0.1.19/bcftools/bcftools
bcftools_1_2=$localToolsDir/samtools/samtools-1.2/bcftools/bcftools-1.2/bcftools
samtools=$localToolsDir/samtools-0.1.19/samtools
samtools_1_2=$localToolsDir/samtools/samtools-1.2/samtools
smmr=`find ${gitDir} -name mmr_for_variant_calling`
varscan2="java -Xmx4g -jar $localToolsDir/varscan2/VarScan.v2.3.7.jar"
somaticSniper=$localToolsDir/somatic-sniper/build/bin/bam-somaticsniper
trimmomatic="java -jar $localToolsDir/Trimmomatic-0.30/trimmomatic-0.30.jar"
mutect="java -Xmx8g -jar $localToolsDir/MuTect-1.1.4/muTect-1.1.4.jar -T MuTect"
GATK="$localToolsDir/jre1.7.0_40/bin/java -jar $localToolsDir/GenomeAnalysisTK-3.3-0/GenomeAna|bowtie_index_tool=<pathTo>/bowtie2-2.1.0/bowtie2-build
lysisTK.jar"
freebayes=$localToolsDir/freebayes/bin/freebayes
jointSNVMix=$localToolsDir/pipeline/bin/jsm.py
pathToJointSNVMix2=$localToolsDir/pipeline/JointSNVMix-0.8-b2/
fasta_generate_regions=$localToolsDir/freebayes/scripts/fasta_generate_regions.py
freebayes_parallel=$localToolsDir/freebayes/scripts/freebayes-parallel
qualimap=<pathTo>/programs/qualimap/qualimap_v2.1.3/qualimap
faToTwoBit=$gitDir/sim_cancer/tools/faToTwoBit
gfServer=$gitDir/sim_cancer/tools/gfServer
gfClient=$gitDir/sim_cancer/tools/gfClient
gemsim_models=$gitDir/sim_cancer/tools/GemSIM_v1.6/models/
fn_model=$gitDir/sim_cancer/tools/GemSIM_v1.6/models/ill100v5_p.gzip
wessim_dir=$gitDir/sim_cancer/tools/Wessim_beta


# reference files
humanRef=<pathTo>/genomes/H_sapiens/hg19.fa
humanRefKaryo=<pathTo>/genomes/H_sapiens/hg19_karyotypic.fa
humanGtf=<pathTo>/genomes/H_sapiens/Homo_sapiens.GRCh37.75.chr.gtf # hg19
fn_genome=$humanRefKaryo
design=S04380219 #:v5 +UTR Nov 2012
bed_file=<pathTo>/agilent/${design}/{$design}_Regions.bed
bed_file_qualimap=${bed_file%.bed}_qualimap.bed
fn_probes=<pathTo>/agilent/$design/${design}_Probes.txt


