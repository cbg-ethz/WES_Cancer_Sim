#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`
source `find ${gitDir} -name paths.sh`
create_cancer_vcfs=`find ${gitDir} -name create_cancer_vcfs.py`
create_ref_from_vcf=`find ${gitDir} -name create_ref_from_vcf`
tools_dir=${gitDir}
source `find ${tools_dir} -name genome.sh`

out_dir=$1 # this folder has to contain the following files: NO_final.vcf, TU_final.vcf, config.tsv
# config.tsv contained this:
# # four levels
# >treeLevelRange
# 0.9 0.5 0.25 0
# >treeStructure
# 1
# 1 1
# 1 1 1 1
# 1 1 1 1 1 1 1 1
# >cygocityT
# 0.5 0.25 0
# >cygocityN
# 0.5 0.25 0


if [ -z "$out_dir" ]; then
	echo "out_dir with NO_final.vcf, TU_final.vcf, and config.tsv not provided"
	exit 0
fi

mkdir -p $out_dir
cd $out_dir


python $create_cancer_vcfs TU_final.vcf NO_final.vcf config.tsv  >> log.log 2>&1

ref=`get_genome human`

for fn_vcf in $out_dir/[0-9]*.vcf
do
	fn_ref_out=${fn_vcf%.vcf}.fa
	fn_map=${fn_vcf%.vcf}_mapping.csv
	command="$create_ref_from_vcf $ref $fn_vcf $fn_ref_out $fn_map >> log.log 2>&1 " 
	echo $command
	eval $command
done

