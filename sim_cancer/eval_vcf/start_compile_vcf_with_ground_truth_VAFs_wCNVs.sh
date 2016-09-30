#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`


# python script with compiles the VAFs for all variants:
python compile_vcf_with_ground_truth_VAFs_wCNVs.py $vcf_sim_dir_wCNVs $base_dir_wCNVs/abundance.txt $base_dir_wCNVs

# to know which variants are in which regions --> important for VAF calculation
bedDir=$vcf_sim_dir_wCNVs 
for subclone_bed in $bedDir/*[01]_TU.bed $bedDir/*[01]_TU_gain.bed ; do
        outVars=${subclone_bed%.bed}_all_variants_regions.bed
        if [ ! -f $outVars ]; then
                command="$intersectBed -f 1.0 -wa -a ~/scratch/arhofman/exome_seq_pipeline_eval/test/all_variants_TEMP.bed -b $subclone_bed > $outVars"
                echo $command
                eval $command
        fi
done

python compile_vcf_with_ground_truth_VAFs_wCNVs.py $vcf_sim_dir_wCNVs $base_dir_wCNVs/abundance.txt $base_dir_wCNVs



