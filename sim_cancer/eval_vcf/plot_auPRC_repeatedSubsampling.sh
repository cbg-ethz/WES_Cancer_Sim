#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`
dir_=`pwd`
namesWoDeep="gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf gatk_SNVs.raw.SOMATIC.vcf  jointSNVMix2_SNVs_Raw.vcf muTect_SNVs_Raw.vcf option__noE_samvar_1_2.SOMATIC.vcf.gz somaticSniper_SNVs_Raw_qual.noComma.vcf bam__varscan2.txt.snp.Somatic_qual.vcf"
cov=50

for eval_dir in $base_dir/alignments_vsn_k20/bam/variants_repeatedSubsampling/eval_10000_160527/; do
	cd $eval_dir
	for fdr in 5 10; do
		# deepSNV
		fn_prc_all=deepSNV.all_auPRC_${fdr}
		if [ ! -f $fn_prc_all ]; then
			for fn_prc in TU.wCont20.final.RG.50perc.*_NO_final.RG.50perc.*_alternative_two.sided.deepSNV.vcf_indel0.eval_auPRC_${fdr}; do
				echo `cat $fn_prc` >> $fn_prc_all 
			done
		fi

		for tool in $namesWoDeep; do
			fn_prc_all=${tool}.all_auPRC_${fdr}
			if [ ! -f $fn_prc_all ]; then
				for fn_prc in TU.wCont20.final.RG.50perc*${tool}*eval_auPRC_${fdr}; do
					echo `cat $fn_prc` >> $fn_prc_all
				done
			fi
		done
	done
done
cd -


for fdr in 5 10; do
	all_auPRC_all_tools="${eval_dir}/deepSNV.all_auPRC_${fdr}  "
	for tool in $namesWoDeep; do
		all_auPRC_all_tools="$all_auPRC_all_tools ${eval_dir}/${tool}.all_auPRC_${fdr}"
	done
	
	command="./plot_auPRC_repeatedSubsampling.R ${eval_dir}/all_auPRC_5_10_repeatedSubsampling_FDR${fdr}.pdf $fdr $all_auPRC_all_tools "
	echo $command
	eval $command
done
	
