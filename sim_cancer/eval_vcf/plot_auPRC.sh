#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`

names="perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf perc.gatk_SNVs.raw.SOMATIC.vcf  perc.jointSNVMix2_SNVs_Raw.vcf perc.muTect_SNVs_Raw.vcf perc.option__noE_samvar_1_2.SOMATIC.vcf.gz perc.somaticSniper_SNVs_Raw_qual.noComma.vcf perc.bam__varscan2.txt.snp.Somatic_qual.vcf"


if true; then
	for eval_dir in $base_dir/alignments_vsn_k20/variants_default/eval_10000_160527/; do
	cd $eval_dir

	for tool in $names 
	do
		for indel in 0
		do
			fn_prc_all=${tool}.all_auPRC_10
			if [ ! -f $fn_prc_all ]; then
				echo > $fn_prc_all
				for cov in 12 25 50 75 100
				do
					if [[ "$tool" == *deep* ]]; then
						 tool="perc_NO_final.RG.${cov}perc_alternative_two.sided.deepSNV.vcf"
					fi 
					fn_prc=TU.wCont20.final.RG.$cov${tool}_indel$indel.eval_auPRC_10 # this is the auPRC when allowing to go down the list until until 10% FDR
					if [ ! -f $fn_prc ]; then 
						continue
					fi
					echo -e "$cov\t`cat $fn_prc`" >> $fn_prc_all
				done
				echo >> $fn_prc_all
				cov=50
				if [[ "$tool" == *deep* ]]; then
					tool="perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf"
				fi
				for cont in 10 20 40 60
				do
					fn_prc=TU.wCont$cont.final.RG.$cov${tool}_indel$indel.eval_auPRC_10 # this is the auPRC when allowing to go down the list until until 10% FDR
					if [ ! -f $fn_prc ]; then 
						continue
					fi
					echo -e "$((cov+$cont-20))\t`cat $fn_prc`" >> $fn_prc_all
				done
			else
				echo "$fn_prc_all exists"
			fi
			fn_prc_all=${tool}.all_auPRC_5
			if [ ! -f $fn_prc_all ]; then
				echo > $fn_prc_all
				for cov in 12 25 50 75 100
				do 
					if [[ "$tool" == *deep* ]]; then
						 tool="perc_NO_final.RG.${cov}perc_alternative_two.sided.deepSNV.vcf"
					fi 
					fn_prc=TU.wCont20.final.RG.$cov${tool}_indel$indel.eval_auPRC_5 # this is the auPRC when allowing to go down the list until until 5% FDR
					if [ ! -f $fn_prc ]; then 
						continue
					fi
					echo -e "$cov\t`cat $fn_prc`" >> $fn_prc_all
				done
				echo >> $fn_prc_all
				cov=50
				if [[ "$tool" == *deep* ]]; then
					tool="perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf"
				fi
				for cont in 10 20 40 60
				do
					fn_prc=TU.wCont$cont.final.RG.$cov${tool}_indel$indel.eval_auPRC_5 # this is the auPRC when allowing to go down the list until until 5% FDR
					if [ ! -f $fn_prc ]; then 
						continue
					fi
					echo -e "$((cov+$cont-20))\t`cat $fn_prc`" >> $fn_prc_all
				done
			else
				echo "$fn_prc_all exists"
			fi
	
		done
	done
	cd -
	

	names="perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf perc.gatk_SNVs.raw.SOMATIC.vcf  perc.jointSNVMix2_SNVs_Raw.vcf perc.muTect_SNVs_Raw.vcf perc.option__noE_samvar_1_2.SOMATIC.vcf.gz perc.somaticSniper_SNVs_Raw_qual.noComma.vcf perc.bam__varscan2.txt.snp.Somatic_qual.vcf"

	namesWodeep="perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf perc.gatk_SNVs.raw.SOMATIC.vcf  perc.jointSNVMix2_SNVs_Raw.vcf perc.muTect_SNVs_Raw.vcf perc.option__noE_samvar_1_2.SOMATIC.vcf.gz perc.somaticSniper_SNVs_Raw_qual.noComma.vcf perc.bam__varscan2.txt.snp.Somatic_qual.vcf"
	
	all_auPRC_all_tools=" ${eval_dir}/perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf.all_auPRC_5 "
	for tool in $namesWodeep
	do
		all_auPRC_all_tools=" $all_auPRC_all_tools ${eval_dir}/${tool}.all_auPRC_5 "
	done
	for tool in $names
	do
		all_auPRC_all_tools=" $all_auPRC_all_tools ${eval_dir}/${tool}.all_auPRC_10 "
	done
	echo $all_auPRC_all_tools
	
	command="${dir_}/plot_auPRC.R ${eval_dir}/all_auPRC_5_10_cov.pdf $all_auPRC_all_tools FALSE"
	echo $command
	eval $command
	
	command="${dir_}/plot_auPRC.R ${eval_dir}/all_auPRC_5_10_cont.pdf $all_auPRC_all_tools TRUE"
	echo $command
	eval $command

	done	
	cd $dir_
fi


## find out which ones are the best performing combinations:
eval_dir=$base_dir/alignments_vsn_k20/variants_combined/eval_10000_160527/
for i in `cat $eval_dir/*auPRC_10 | sort -n | tail -n 10`; do file=`grep -l $i $eval_dir/*auPRC_10`; echo -e "`basename $file`\t$i"; done
echo -e "\n"
for i in `cat $eval_dir/*auPRC_5 | sort -n | tail -n 10`; do file=`grep -l $i $eval_dir/*auPRC_5`; echo -e "`basename $file`\t$i"; done
## 
## To get the top 10 of auPRC5 and auPRC10 in the right format for the latex table
for i in `cat $eval_dir/*auPRC_10 | sort -n | tail -n 10`; do file=`grep -l $i $eval_dir/*auPRC_10`; tools=`basename $file | sed 's/combinedlistbetter_//' | sed 's/_prod.vcf_indel0.eval_auPRC_10//' | tr '_' ',' | sed 's/varscan2/VarScan2/' | sed 's/gatk,SNVs/GATK UG/' | sed 's/muTect/MuTect/' | sed 's/jointSNVMix2/JointSNVMix2/' | sed 's/gatkHPCaller/GATK HP/' | sed 's/samvar,1,2/SAMtools/'`; file5perc=${file%auPRC_10}auPRC_5; auprc5=`cat $file5perc`; echo -e "$tools\t& $auprc5 &\t$i "; done | tac

for i in `cat $eval_dir/*auPRC_5 | sort -n | tail -n 10`; do file=`grep -l $i $eval_dir/*auPRC_5`; tools=`basename $file | sed 's/combinedlistbetter_//' | sed 's/_prod.vcf_indel0.eval_auPRC_5//' | tr '_' ','  | sed 's/varscan2/VarScan2/' | sed 's/gatk,SNVs/GATK UG/' | sed 's/muTect/MuTect/' | sed 's/jointSNVMix2/JointSNVMix2/' | sed 's/gatkHPCaller/GATK HP/' | sed 's/samvar,1,2/SAMtools/'`; file10perc=${file%auPRC_5}auPRC_10; auprc10=`cat $file10perc`; echo -e "$tools\t& $i &\t$auprc10 "; done | tac
##



