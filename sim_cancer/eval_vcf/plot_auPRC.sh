#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`

names="perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf perc.gatk_SNVs.raw.SOMATIC.vcf  perc.jointSNVMix2_SNVs_Raw.vcf perc.muTect_SNVs_Raw.vcf perc.option__noE_samvar_1_2.SOMATIC.vcf.gz perc.sinvict_TU.wCont.final.RG.perc_vs_NO_final.RG.perc_somatic.vcf  perc.somaticSniper_SNVs_Raw_qual.noComma.vcf perc.bam__varscan2.txt.snp.Somatic_qual.vcf"


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
					if [[ "$tool" == *sinvict* ]]; then
						fn_prc=sinvict_TU.wCont20.final.RG.${cov}perc_vs_NO_final.RG.${cov}perc_somatic.vcf_indel${indel}.eval_auPRC_10
					fi
					if [ ! -f $fn_prc ]; then 
						echo "$fn_prc does not exist"
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
					if [[ "$tool" == *sinvict* ]]; then
						fn_prc=sinvict_TU.wCont${cont}.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel${indel}.eval_auPRC_10
					fi
					if [ ! -f $fn_prc ]; then 
						echo "$fn_prc does not exist"
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
					if [[ "$tool" == *sinvict* ]]; then
						fn_prc=sinvict_TU.wCont20.final.RG.${cov}perc_vs_NO_final.RG.${cov}perc_somatic.vcf_indel${indel}.eval_auPRC_5
					fi
					if [ ! -f $fn_prc ]; then 
						echo "$fn_prc does not exist"
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
					if [[ "$tool" == *sinvict* ]]; then
						fn_prc=sinvict_TU.wCont${cont}.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel${indel}.eval_auPRC_5
					fi
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
	
	names="perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf perc.gatk_SNVs.raw.SOMATIC.vcf  perc.jointSNVMix2_SNVs_Raw.vcf perc.muTect_SNVs_Raw.vcf perc.option__noE_samvar_1_2.SOMATIC.vcf.gz perc.sinvict_TU.wCont.final.RG.perc_vs_NO_final.RG.perc_somatic.vcf  perc.somaticSniper_SNVs_Raw_qual.noComma.vcf perc.bam__varscan2.txt.snp.Somatic_qual.vcf"

	namesWodeep="perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf perc.gatk_SNVs.raw.SOMATIC.vcf  perc.jointSNVMix2_SNVs_Raw.vcf perc.muTect_SNVs_Raw.vcf perc.option__noE_samvar_1_2.SOMATIC.vcf.gz perc.sinvict_TU.wCont.final.RG.perc_vs_NO_final.RG.perc_somatic.vcf  perc.somaticSniper_SNVs_Raw_qual.noComma.vcf perc.bam__varscan2.txt.snp.Somatic_qual.vcf"
	
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
	
	command="${dir_}/sim_cancer/eval_vcf/plot_auPRC.R ${eval_dir}/all_auPRC_5_10_cov.pdf $all_auPRC_all_tools FALSE"
	echo $command
	eval $command
	
	command="${dir_}/sim_cancer/eval_vcf/plot_auPRC.R ${eval_dir}/all_auPRC_5_10_cont.pdf $all_auPRC_all_tools TRUE"
	echo $command
	eval $command

	done	
	cd $dir_
fi
exit 0

## find out which ones are the best performing combinations:
eval_dir=$base_dir/alignments_vsn_k20/variants_combined/eval_10000_160527/
eval_dir_sinvict=$base_dir/alignments_vsn_k20/variants_combined_sinvict/eval_10000_160527/
for i in `cat $eval_dir/*auPRC_10 | sort -n | tail -n 10`; do file=`grep -l $i $eval_dir/*auPRC_10`; echo -e "`basename $file`\t$i"; done
echo -e "\n"
for i in `cat $eval_dir/*auPRC_5 | sort -n | tail -n 10`; do file=`grep -l $i $eval_dir/*auPRC_5`; echo -e "`basename $file`\t$i"; done
## 
## To get the top 10 of auPRC5 and auPRC10 in the right format for the latex table
for i in `cat $eval_dir/*auPRC_10 | sort -n | tail -n 10`; do file=`grep -l $i $eval_dir/*auPRC_10`; tools=`basename $file | sed 's/combinedlistbetter_[456]_//' | sed 's/_prod.vcf_indel0.eval_auPRC_10//' | tr '_' ',' | sed 's/varscan2/VarScan2/' | sed 's/gatk,SNVs/GATK UG/' | sed 's/muTect/MuTect/' | sed 's/jointSNVMix2/JointSNVMix2/' | sed 's/gatkHPCaller/GATK HP/' | sed 's/samvar,1,2/SAMtools/' | sed 's/sinvict/SiNVICT/'`; file5perc=${file%auPRC_10}auPRC_5; auprc5=`cat $file5perc`; echo -e "$tools\t& num{$auprc5} &num{$i} "; done | tac

for i in `cat $eval_dir_sinvict/*auPRC_10 | sort -n | tail -n 11`; do file=`grep -l $i $eval_dir_sinvict/*auPRC_10 | awk 'NR==1'`; tools=`basename $file | sed 's/combinedlistbetter_[456]_//' | sed 's/_prod.vcf_indel0.eval_auPRC_10//' | tr '_' ',' | sed 's/deep/deepSNV/' | sed 's/varsc/VarScan2/' | sed 's/UG/GATK UG/' | sed 's/mutec/MuTect/' | sed 's/joint/JointSNVMix2/' | sed 's/HP/GATK HP/' | sed 's/samv,1,2/SAMtools/' | sed 's/sinv/SiNVICT/' | sed 's/snipe/somaticSniper/'`; file5perc=${file%auPRC_10}auPRC_5; auprc5=`cat $file5perc`; filename=`basename $file`; if [[ "$filename" == *sinv* || "$filename" == *all* ]]; then echo -e "$tools\t&num{$auprc5}&num{$i}"; fi; done | tac



for i in `cat $eval_dir/*auPRC_5 | sort -n | tail -n 10`; do file=`grep -l $i $eval_dir/*auPRC_5`; tools=`basename $file | sed 's/combinedlistbetter_//' | sed 's/_prod.vcf_indel0.eval_auPRC_5//' | tr '_' ','  | sed 's/varscan2/VarScan2/' | sed 's/gatk,SNVs/GATK UG/' | sed 's/muTect/MuTect/' | sed 's/jointSNVMix2/JointSNVMix2/' | sed 's/gatkHPCaller/GATK HP/' | sed 's/samvar,1,2/SAMtools/' | sed 's/sinvict/SiNVICT/'`; file10perc=${file%auPRC_5}auPRC_10; auprc10=`cat $file10perc`; echo -e "$tools\t& $i &\t$auprc10 "; done | tac

for i in `cat $eval_dir_sinvict/*auPRC_5 | sort -n | tail -n 11`; do file=`grep -l $i $eval_dir_sinvict/*auPRC_5 | awk 'NR==1'`; tools=`basename $file | sed 's/combinedlistbetter_[456]_//' | sed 's/_prod.vcf_indel0.eval_auPRC_10//' | tr '_' ',' | sed 's/deep/deepSNV/' | sed 's/varsc/VarScan2/' | sed 's/UG/GATK UG/' | sed 's/mutec/MuTect/' | sed 's/joint/JointSNVMix2/' | sed 's/HP/GATK HP/' | sed 's/samv,1,2/SAMtools/' | sed 's/sinv/SiNVICT/' | sed 's/snipe/somaticSniper/'`; file10perc=${file%auPRC_5}auPRC_10; auprc10=`cat $file10perc`; filename=`basename $file`; if [[ "$filename" == *sinv* || "$filename" == *all* ]]; then echo -e "$tools\t&num{i}&num{$auprc10}"; fi; done | tac
##



