#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`


echo "${currScriptDir}/../simulate_WES_cancerSample/repeated_subsampling.sh $base_dir/alignments_vsn_k20/bam/"
echo "${currScriptDir}/run_variant_callers.sh $base_dir/alignments_vsn_k20/bam/repeatedSubsampling/ $base_dir/alignments_vsn_k20/bam/  perc.[0-9]\*.bam"

for out_folder in `ls -lh $base_dir/alignments*/bam/*.bam | awk '{print $9}'  | sed 's!/bam/! !' | awk '{print $1}' | uniq` 
do
	echo -e "\n"
	echo "${currScriptDir}/../simulate_WES_cancerSample/process_alignments.sh $out_folder/bam true" 
	echo "${currScriptDir}/../simulate_WES_cancerSample/process_true_alignments.sh $out_folder/bam" 
done

for out_folder in `ls -lh $base_dir/alignments*/bam/*.bam | awk '{print $9}'  | sed 's!/bam/! !' | awk '{print $1}' | uniq`
do
	echo -e "\n"
	echo "${currScriptDir}/run_variant_callers.sh $out_folder/bam $out_folder " 
	echo "${currScriptDir}/run_variant_callers.sh $out_folder/bam $out_folder perc.true.bam" 
	echo "loc_realign=true"
	echo "${currScriptDir}/run_variant_callers.sh $out_folder/bam $out_folder"
	echo "${currScriptDir}/run_variant_callers.sh $out_folder/bam $out_folder perc.true.bam"

	if [[ $out_folder == *vsn_k20* ]]; then
		echo "${currScriptDir}/tune_variant_callers.sh $out_folder/bam $out_folder"
		echo "${currScriptDir}/tune_variant_callers.sh $out_folder/bam $out_folder perc.true.bam"
	fi
done

echo "${currScriptDir}/../simulate_WES_cancerSample/check_coverage.sh $base_dir/alignments_vsn_k20/bam $base_dir/alignments_vsn_k20"

