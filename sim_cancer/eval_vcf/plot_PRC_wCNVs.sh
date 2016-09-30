#!/bin/bash


currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`

path2PR="$base_dir_wCNVs/alignments_vsn_k20/variants_default/eval_10000_160527/"

declare -a tools=('_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf' '.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf' '.gatk_SNVs.raw.SOMATIC.vcf' '.jointSNVMix2_SNVs_Raw.vcf' '.muTect_SNVs_Raw.vcf' '.option__noE_samvar_1_2.SOMATIC.vcf.gz' 'sinvict' '.somaticSniper_SNVs_Raw_qual.noComma.vcf' '.bam__varscan2.txt.snp.Somatic_qual.vcf')
declare -a titles=('deepSNV' 'GATK_HP' 'GATK_UG' 'JointSNVMix2' 'MuTect' 'SAMtools' 'SiNVICT' 'SomaticSniper' 'VarScan2')


beginnig="TU.wCont20.final.RG.50perc"
ending="_indel0.PRcurve"


allFiles=""
for fileIDX in `seq 1 ${#tools[@]}`; do
	myPRcurveFile=$path2PR/${beginnig}${tools[fileIDX-1]}${ending}
	if [[ "${tools[fileIDX-1]}" == "sinvict"  ]]; then
		myPRcurveFile=$path2PR/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.PRcurve
	fi
	allFiles="$allFiles $myPRcurveFile "
	command="./plot_PRC.R $path2PR $myPRcurveFile ${titles[fileIDX-1]}"
	echo $command
	eval $command
done

./plot_PRC_altogether.R $allFiles



