#!/bin/bash


currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`

### filter only somatic mutations #locally
varDir=$base_dir/alignments_vsn_k20/variants_default/
var_dir=$varDir
cont=20
perc=50
true_tag=""

outDir=$varDir/binomTest_160527
mkdir -p $outDir

for mut_list in $var_dir/TU.wCont${cont}.final.RG.${perc}perc${true_tag}.jointSNVMix2_SNVs_Raw.vcf
do
	echo $mut_list
	fileName=`basename ${mut_list%.gz}`
	outFile=$outDir/${fileName}_cpp
	echo "JointSNVMix2"
	$currScriptDir/filter_somatic_binom_test_JointSNVMix2 $mut_list ${outFile}
done


