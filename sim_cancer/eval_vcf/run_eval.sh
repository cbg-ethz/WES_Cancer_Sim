#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`

wCNVs=$1

if [ "$1" == "wCNVs" ]; then
	out_dir=$base_dir_wCNVs
	vcf_dir=$vcf_sim_dir_wCNVs
else
	out_dir=$base_dir
	vcf_dir=$vcf_sim_dir
fi

echo $out_dir
echo $vcf_dir

mkdir -p $out_dir

#####
## First:
# - transform all the varscan and somatic sniper VCF's into the *_qual.vcf's with the quality of a variant in the sixth column -> Only then, a meaningful evaluation can be done
# - filter out all germline mutations from GATK UG, HAPCaller, SAMtools. These tools also report germline mutations and give them a high score -> this obscurs the comparison
# 	makes them look worse as deepSNV, JSM2, VarScan2, SomaticSniper, who explicitly sort by some sort of somatic score. -> We run them separate on TU and NO. We divide the TU vcf up into germline and somatic, based on what can be found also in the NO vcf. If there are several alternate alleles, we only write out those that are not also in the normal sample!
# - format GATKHPCaller: sometimes the FORMAT field of GATKHPCaller is malformed -> the vcf parser crashes. Therefore we reformat it.
# - for somaticSniper, we do the filter script filterOutGermlineComma_somSniper.sh because: There are cases like this A -> C,G, where C is a germline mutation, and G is the actual somatic mutation. In the evaluation script, these would both be understood as somatic mutation --> false positive because germline. In the filter script, we take the genotype from TU and NO, and only output those alternate alleles that represent somatic mutations
#####

max_num_jobs=50
#eval_file_ending=FN_absolute_counts
eval_file_ending=SN

queue="mpi01.q"
cnt_jobs_TODO=0
cnt_jobs_total=0
source `find ${gitDir} -name submit.sh` 
source `find ${gitDir} -name genome.sh`
max_cnt_normal=10000
evaluate_vcf_tool_path=$currScriptDir

wait_for_jobs eval_SNV --max-num $max_num_jobs

# # ##############################################
# # ## Assessing the "PASS" label of mutect:
# # ## cd $out_dir
# for i in  ./alignments_vsn_k20/variants_default/TU.wCont20.final.RG.*perc.muTect_SNVs_Raw.vcf; do echo $i; if [ ! -f ${i%_Raw.vcf}.noReject.vcf ]; then cat $i | awk '$7!="REJECT"' > ${i%_Raw.vcf}.noReject.vcf;  fi; done
# # ## for i in ./alignments_vsn_k20/variants_default/TU.wCont20.final.RG.*perc.muTect_SNVs_Raw.vcf; do echo $i; if [ ! -f ${i%_Raw.vcf}.noPass.vcf ]; then cat $i | awk '$7!="PASS"' > ${i%_Raw.vcf}.noPass.vcf;  fi; done
# # ## cd -
# var_dir=$out_dir/alignments_vsn_k20/variants_default/
# true_tag=""
# for perc in 12 25 50 75 100; do
# 	cont=20
# 	fn_bam=$var_dir/../bam/TU.wCont${cont}.final.RG.${perc}perc.bam
# 	fn_bam_normal=$var_dir/../bam/NO_final.RG.${perc}perc.bam 
# 	mutect="$var_dir/TU.wCont$cont*${perc}perc${true_tag}.muTect_SNVs.noReject.vcf"
# 	#mutect="$var_dir/TU.wCont$cont*${perc}perc${true_tag}.muTect_SNVs.noPass.vcf"
# 	tool=$mutect
# 	eval_dir=$var_dir/eval_${noY_}${max_cnt_normal}_160527
# 	fn_out=$eval_dir/logs/`basename $tool`_${RANDOM}.o
# 	mkdir -p $eval_dir
# 	fn_SN=$eval_dir/`basename $tool`_indel0.eval_${eval_file_ending}
# 	command="$evaluate_vcf_tool_path/evaluate_vcf --out-dir $out_dir --vcf-dir $vcf_dir --eval-dir $eval_dir --max-cnt-normal $max_cnt_normal $tool --bam $fn_bam --bam-normal $fn_bam_normal --only-sn-auprc --pathToBam $out_dir --pathToBed $pathToBed"
# 	submit "$command" --tag eval_SNV --log-dir $eval_dir/logs --queue $queue --fn-out $fn_out
# 	echo "$command"
# done
# exit 0
# # ##############################################

# ##################################################
# ## Assessing the different levels of SiNVICT
# ## cat ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf | grep -v ^# | awk '{print $6}' | sort | uniq -c
# ## cat ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf | grep ^# > ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_HEADER
# ## cat ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf | grep -v ^# | awk '$6 == 4' > ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.Level4.vcf_unheadered
# ## cat ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf | grep -v ^# | awk '$6 == 2' > ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.Level2.vcf_unheadered
# ## cat ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf | grep -v ^# | awk '$6 == 1' > ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.Level1.vcf_unheadered
# ## cat ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_HEADER ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.Level4.vcf_unheadered > ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.Level4.vcf
# ## cat ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_HEADER ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.Level2.vcf_unheadered > ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.Level2.vcf
# ## cat ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_HEADER ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.Level1.vcf_unheadered > ./alignments_vsn_k20/variants_default/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.Level1.vcf
# var_dir=$out_dir/alignments_vsn_k20/variants_default/
# true_tag=""
# perc=50
# cont=20
# fn_bam=$var_dir/../bam/TU.wCont${cont}.final.RG.${perc}perc.bam
# fn_bam_normal=$var_dir/../bam/NO_final.RG.${perc}perc.bam
# for level in 1 2 4; do
# 	tool=$var_dir/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.Level${level}.vcf
# 	eval_dir=$var_dir/eval_${noY_}${max_cnt_normal}_160527	
# 	fn_out=$eval_dir/logs/`basename $tool`_${RANDOM}.o	
# 	mkdir -p $eval_dir
# 	command="$evaluate_vcf_tool_path/evaluate_vcf --out-dir $out_dir --vcf-dir $vcf_dir --eval-dir $eval_dir --max-cnt-normal $max_cnt_normal $tool --bam $fn_bam --bam-normal $fn_bam_normal --only-sn-auprc --pathToBam $out_dir --pathToBed $pathToBed"
# 	submit "$command" --tag eval_SNV --log-dir $eval_dir/logs --queue $queue --fn-out $fn_out
# 	echo "$command"
# done
# exit 0
# ##################################################


## Set the runs you wish to evaluate to 'true'
default=false
tuned=false
locRe=false
combis=false
overlps=true
binomT=false
repeatS=false

if $default; then
	# first the default runs:
	for align_version in vsn_k20 sn_k1 vsn_k20 
	do
	var_dir=$out_dir/alignments_${align_version}/variants_default/
	for true_tag in "" ; do

		for perc in 50 12 25 75 100; do
			for cont in 20 10 40 60; do
				if [ $cont != 20 -a $perc != 50 ]; then
					continue
				fi

				if [[ $out_dir == *Aneuploidy*  || $align_version == "sn_k1" ]]; then
			                if [ ! $perc == 50 -o ! $cont == 20  ]; then
			                        continue
			                fi
			        fi

				fn_bam=$var_dir/../bam/TU.wCont${cont}.final.RG.${perc}perc.bam
				fn_bam_normal=$var_dir/../bam/NO_final.RG.${perc}perc.bam 

				jointSNV="$var_dir/TU.wCont${cont}.final.RG.${perc}perc${true_tag}.jointSNVMix2_SNVs_Raw.vcf"
				varscan="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.bam__varscan2.txt.snp.Somatic_qual.vcf"
				sniper="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.somaticSniper_SNVs_Raw_qual.noComma.vcf"
				sam_1_2="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.option__noE_samvar_1_2.SOMATIC.vcf.gz"
				gatk="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.gatk_SNVs.raw.SOMATIC.vcf"
				gatkHPCaller="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf"
				deepSNV="$var_dir/TU.wCont${cont}.final.RG.${perc}perc${true_tag}_NO_final.RG.${perc}perc${true_tag}_alternative_two.sided.deepSNV.vcf"
 				mutect="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.muTect_SNVs_Raw.vcf"
				sinvict="$var_dir/sinvict_TU.wCont${cont}.final.RG.${perc}perc_vs_NO_final.RG.${perc}perc_somatic.vcf"
			
				for tool in $sniper $jointSNV $sam_1_2 $varscan $gatk $deepSNV $mutect $gatkHPCaller $sinvict; do
					vcf_name=`basename $tool`

					cnt_jobs_total=$((cnt_jobs_total+1))
	
					eval_dir=`dirname $tool`/eval_${noY_}${max_cnt_normal}_160527
					mkdir -p $eval_dir
					mkdir -p $eval_dir/logs
					fn_out=$eval_dir/logs/`basename $tool`_${RANDOM}.o

					if true; then
						fn_SN=$eval_dir/`basename $tool`_indel0.eval_${eval_file_ending}

						exec_creation_time=`stat -c %Y $evaluate_vcf_tool_path/evaluate_vcf `
						fn_SN_creation_time=`stat -c %Y $fn_SN`	
						if [ -f $fn_SN ]; then 
							fn_SN_creation_time=`stat -c %Y $fn_SN`	
							if [ $fn_SN_creation_time -lt $exec_creation_time ]; then
								echo $fn_SN exists, and is older
							else
								echo $fn_SN exists, and is newer
								continue
							fi
						else
							echo $fn_SN does not exist
						fi
					fi

					if [ $cont == 20 -a $perc == 50 ]; then 
						command="$evaluate_vcf_tool_path/evaluate_vcf --out-dir $out_dir --vcf-dir $vcf_dir --eval-dir $eval_dir --max-cnt-normal $max_cnt_normal $tool --bam $fn_bam --bam-normal $fn_bam_normal --pathToBam $out_dir --pathToBed $pathToBed"
					else
						command="$evaluate_vcf_tool_path/evaluate_vcf --out-dir $out_dir --vcf-dir $vcf_dir --eval-dir $eval_dir --max-cnt-normal $max_cnt_normal $tool --bam $fn_bam --bam-normal $fn_bam_normal --only-sn-auprc --pathToBam $out_dir --pathToBed $pathToBed " 
					fi
					echo "$command"
					submit "$command" --tag eval_SNV --log-dir $eval_dir/logs --queue $queue --fn-out $fn_out
					cnt_jobs_TODO=$((cnt_jobs_TODO+1))
					wait_for_jobs eval_SNV --max-num $max_num_jobs
				done
			done
		done
		wait_for_jobs eval_SNV --max-num $max_num_jobs
	done
	done
	echo "cnt_jobs_TODO=$cnt_jobs_TODO"
fi


if $tuned; then	
	# now the tuned runs:
	for var_dir in $out_dir/alignments_vsn_k20/variants_tuned/; do
		perc=50
		cont=20
		fn_bam=$var_dir/../bam/TU.wCont${cont}.final.RG.${perc}perc.bam
		fn_bam_normal=$var_dir/../bam/NO_final.RG.${perc}perc.bam 
	
		for vcf in `ls $var_dir/*joint*.vcf | grep -v combined` `ls $var_dir/*deep*.vcf | grep -v combined` `ls $var_dir/*samvar*SOMATIC*.gz` `ls $var_dir/*varscan*qual.vcf | grep -v combined` `ls $var_dir/sinvict*.vcf | grep -v combined` ; do
			cnt_jobs_total=$((cnt_jobs_total+1))
		
			eval_dir=$var_dir/eval_${noY_}${max_cnt_normal}_160527
			mkdir -p $eval_dir
	
			if true; then
				fn_SN=$eval_dir/`basename $vcf`_indel0.eval_${eval_file_ending}
	
				exec_creation_time=`stat -c %Y $evaluate_vcf_tool_path/evaluate_vcf `
				fn_SN_creation_time=`stat -c %Y $fn_SN`	
				if [ -f $fn_SN ]; then 
					fn_SN_creation_time=`stat -c %Y $fn_SN`	
					if [ $fn_SN_creation_time -lt $exec_creation_time ]; then
						echo $fn_SN exists, and is older
					else
						echo $fn_SN exists, and is newer
						continue
					fi
				else
					echo $fn_SN does not exist
				fi
			fi
	
			command="$evaluate_vcf_tool_path/evaluate_vcf --out-dir $out_dir --vcf-dir $vcf_dir --eval-dir $eval_dir --max-cnt-normal $max_cnt_normal $vcf --bam $fn_bam --bam-normal $fn_bam_normal  --only-sn-auprc --pathToBam $out_dir --pathToBed $pathToBed "
			echo "$command"
			submit "$command" --tag eval_SNV --log-dir $eval_dir/logs --queue $queue
			cnt_jobs_TODO=$((cnt_jobs_TODO+1))
			wait_for_jobs eval_SNV --max-num $max_num_jobs
		done
	done
	echo "cnt_jobs_TODO=$cnt_jobs_TODO"
#exit 0
fi


if $locRe; then
	# now the local realignment runs
	for var_dir in $out_dir/alignments_vsn_k20/variants_default/; do
		perc=50
		cont=20
		fn_bam=$var_dir/../bam/TU.wCont${cont}.final.RG.${perc}perc.bam
		fn_bam_normal=$var_dir/../bam/NO_final.RG.${perc}perc.bam 
	
		for true_tag in ""; do
			jointSNV="$var_dir/TU.wCont${cont}.final.RG.${perc}perc${true_tag}.locRealign.jointSNVMix2_SNVs_Raw.vcf"
			varscan="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.locRealign.bam__varscan2.txt.snp.Somatic_qual.vcf"
			sniper="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.locRealign.somaticSniper_SNVs_Raw_qual.noComma.vcf"
			sam_1_2="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.locRealign.option__noE_samvar_1_2.SOMATIC.vcf.gz"
			gatk="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.locRealign.gatk_SNVs.raw.SOMATIC.vcf"
			gatkHPCaller="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.locRealign.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf"
	 		deepSNV="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.locRealign_NO_final.RG.${perc}perc.locRealign_alternative_two.sided.deepSNV.vcf"
	 		mutect="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.locRealign.muTect_SNVs_Raw.vcf"
			sinvict="$var_dir/sinvict_TU.wCont${cont}.final.RG.${perc}perc.locRealign_vs_NO_final.RG.${perc}perc.locRealign_somatic.vcf"
	
	
			for vcf in $jointSNV $varscan $sniper $sam_1_2 $gatk $gatkHPCaller $deepSNV $mutect $sinvict ; do
				cnt_jobs_total=$((cnt_jobs_total+1))
		
				eval_dir=$var_dir/eval_${noY_}${max_cnt_normal}_160527
				mkdir -p $eval_dir
	
				if true; then
					fn_SN=$eval_dir/`basename $vcf`_indel0.eval_${eval_file_ending}
	
					exec_creation_time=`stat -c %Y $evaluate_vcf_tool_path/evaluate_vcf `
					fn_SN_creation_time=`stat -c %Y $fn_SN`	
					if [ -f $fn_SN ]; then 
						fn_SN_creation_time=`stat -c %Y $fn_SN`	
						if [ $fn_SN_creation_time -lt $exec_creation_time ]; then
							echo $fn_SN exists, and is older
						else
							echo $fn_SN exists, and is newer
							continue
						fi
					else
						echo $fn_SN does not exist
					fi
				fi
	
				command="$evaluate_vcf_tool_path/evaluate_vcf --out-dir $out_dir --vcf-dir $vcf_dir --eval-dir $eval_dir --max-cnt-normal $max_cnt_normal $vcf --bam $fn_bam --bam-normal $fn_bam_normal  --only-sn-auprc --pathToBam $out_dir --pathToBed $pathToBed "
				echo "$command"
				submit "$command" --tag eval_SNV --log-dir $eval_dir/logs --queue $queue
				cnt_jobs_TODO=$((cnt_jobs_TODO+1))
				wait_for_jobs eval_SNV --max-num $max_num_jobs
			done
		done
	done
fi

	
if $combis; then	
	# now the combinations
	#for curr_dir in $out_dir/alignments_vsn_k20/variants_combined/; do
	for curr_dir in $out_dir/alignments_vsn_k20/variants_combined_sinvict/; do
		perc=50
		cont=20
		fn_bam=$curr_dir/../bam/TU.wCont${cont}.final.RG.${perc}perc.bam
		fn_bam_normal=$curr_dir/../bam/NO_final.RG.${perc}perc.bam
	
		for vcf in `find $curr_dir -name combinedlistbetter_*prod.vcf | grep -v eval_1000 | grep -v oldComb`; do
	
			echo $vcf
	
			cnt_jobs_total=$((cnt_jobs_total+1))
		
			eval_dir=`dirname $vcf`/eval_${noY_}${max_cnt_normal}_160527
			mkdir -p $eval_dir
			fn_out=$eval_dir/logs/`basename $vcf`_${RANDOM}.o
			if [ -f $fn_out ]; then
				continue
			fi
			if true; then
				fn_SN=$eval_dir/`basename $vcf`_indel0.eval_${eval_file_ending}
	
				exec_creation_time=`stat -c %Y $evaluate_vcf_tool_path/evaluate_vcf `
				fn_SN_creation_time=`stat -c %Y $fn_SN`	
				if [ -f $fn_SN ]; then 
					fn_SN_creation_time=`stat -c %Y $fn_SN`	
					if [ $fn_SN_creation_time -lt $exec_creation_time ]; then
						echo $fn_SN exists, and is older
					else
						echo $fn_SN exists, and is newer
						continue
					fi
				else
					echo $fn_SN does not exist
				fi
			fi
			command="$evaluate_vcf_tool_path/evaluate_vcf --out-dir $out_dir --vcf-dir $vcf_dir --eval-dir $eval_dir --max-cnt-normal $max_cnt_normal $vcf --bam $fn_bam --bam-normal $fn_bam_normal  --only-sn-auprc --pathToBam $out_dir --pathToBed $pathToBed"
			echo "$command"
			submit "$command" --tag eval_SNV --log-dir $eval_dir/logs --queue $queue --fn-out $fn_out
			cnt_jobs_TODO=$((cnt_jobs_TODO+1))
			wait_for_jobs eval_SNV --max-num $max_num_jobs
		done
	done
#exit 0
fi

if $overlps; then
	## now the overlaps
	#for var_dir in $overlapDir/specific_overlaps/overlaps; do
	#for var_dir in $overlapDir/SNV_lists_fdr01_1percent/overlaps/allOverlapsVennTop5ToolsExtractLists/ $overlapDir/SNV_lists_fdr05_5percent/overlaps/allOverlapsVennTop5ToolsExtractLists/ $overlapDir/SNV_lists_fdr10_10percent/overlaps/allOverlapsVennTop5ToolsExtractLists/ ; do
	for var_dir in $overlapDir/all_pairwise_overlaps/pairwise_overlaps/; do

		algnmt_vsn_k20=$out_dir/alignments_vsn_k20/
		perc=50
		cont=20
		fn_bam=$algnmt_vsn_k20/bam/TU.wCont${cont}.final.RG.${perc}perc.bam
		fn_bam_normal=$algnmt_vsn_k20/bam/NO_final.RG.${perc}perc.bam 

		for vcf in $var_dir/*original*.vcf; do
			echo $vcf

			cnt_jobs_total=$((cnt_jobs_total+1))
		
			eval_dir=$var_dir/eval_${noY_}${max_cnt_normal}_160527
			fn_out=$eval_dir/logs/`basename $vcf`_${RANDOM}.o
			
			mkdir -p $eval_dir
	
			if true; then
				fn_SN=$eval_dir/`basename $vcf`_indel0.eval_${eval_file_ending}
	
				exec_creation_time=`stat -c %Y $evaluate_vcf_tool_path/evaluate_vcf `
				fn_SN_creation_time=`stat -c %Y $fn_SN`	
				if [ -f $fn_SN ]; then 
					fn_SN_creation_time=`stat -c %Y $fn_SN`	
					if [ $fn_SN_creation_time -lt $exec_creation_time ]; then
						echo $fn_SN exists, and is older
					else
						echo $fn_SN exists, and is newer
						continue
					fi
				else
					echo $fn_SN does not exist
				fi
			fi
			command="$evaluate_vcf_tool_path/evaluate_vcf --out-dir $out_dir --vcf-dir $vcf_dir --eval-dir $eval_dir --max-cnt-normal $max_cnt_normal $vcf --bam $fn_bam --bam-normal $fn_bam_normal  --only-sn-auprc --pathToBam $out_dir --pathToBed $pathToBed "
			echo "$command"
			submit "$command" --tag eval_SNV --log-dir $eval_dir/logs --queue $queue --fn-out $fn_out
			cnt_jobs_TODO=$((cnt_jobs_TODO+1))
			wait_for_jobs eval_SNV --max-num $max_num_jobs
		done
	done
#exit 0
fi


if $binomT; then
	## Now the binomial test ones
	for var_dir in $out_dir/alignments_vsn_k20/variants_default/; do
		perc=50
		cont=20
		true_tag=""
		fn_bam=$var_dir/../bam/TU.wCont${cont}.final.RG.${perc}perc.bam
		fn_bam_normal=$var_dir/../bam/NO_final.RG.${perc}perc.bam 
		
		jointSNV="$var_dir/TU.wCont${cont}.final.RG.${perc}perc${true_tag}.jointSNVMix2_SNVs_Raw.vcf"
		varscan="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.bam__varscan2.txt.snp.Somatic_qual.vcf"
		sniper="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.somaticSniper_SNVs_Raw_qual.noComma.vcf"
		sam_1_2="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.option__noE_samvar_1_2.SOMATIC.vcf.gz"
		gatk="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.gatk_SNVs.raw.SOMATIC.vcf"
		gatkHPCaller="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf"
		deepSNV="$var_dir/TU.wCont${cont}.final.RG.${perc}perc${true_tag}_NO_final.RG.${perc}perc${true_tag}_alternative_two.sided.deepSNV.vcf"
 		mutect="$var_dir/TU.wCont$cont.final.RG.${perc}perc${true_tag}.muTect_SNVs_Raw.vcf"
		sinvict="$var_dir/sinvict_TU.wCont${cont}.final.RG.${perc}perc_vs_NO_final.RG.${perc}perc_somatic.vcf"
			
		for tool in $sniper $jointSNV $sam_1_2 $varscan $gatk $deepSNV $mutect $gatkHPCaller $sinvict; do
			vcf_name=`basename $tool`
			vcf=$tool
			cnt_jobs_total=$((cnt_jobs_total+1))
		
			eval_dir=$var_dir/eval_${noY_}${max_cnt_normal}_binomTest_160527
			fn_out=$eval_dir/logs/`basename $vcf`_${RANDOM}.o
			mkdir -p $eval_dir
	
			if true; then
				fn_SN=$eval_dir/`basename $vcf`_indel0.eval_${eval_file_ending}
	
				exec_creation_time=`stat -c %Y $evaluate_vcf_tool_path/evaluate_vcf `
				fn_SN_creation_time=`stat -c %Y $fn_SN`	
				if [ -f $fn_SN ]; then 
					fn_SN_creation_time=`stat -c %Y $fn_SN`	
					if [ $fn_SN_creation_time -lt $exec_creation_time ]; then
						echo $fn_SN exists, and is older
					else
						echo $fn_SN exists, and is newer
						continue
					fi
				else
					echo $fn_SN does not exist
				fi
			fi
	
			command="$evaluate_vcf_tool_path/evaluate_vcf --out-dir $out_dir --vcf-dir $vcf_dir --eval-dir $eval_dir --max-cnt-normal $max_cnt_normal $tool --bam $fn_bam --bam-normal $fn_bam_normal --only-sn-auprc --binom-test --pathToBam $out_dir --pathToBed $pathToBed "
			echo "$command"
			submit "$command" --tag eval_SNV --log-dir $eval_dir/logs --fn-out $fn_out --queue $queue
			cnt_jobs_TODO=$((cnt_jobs_TODO+1))
			wait_for_jobs eval_SNV --max-num $max_num_jobs
		done
	done
	echo "cnt_jobs_TODO=$cnt_jobs_TODO"

	## from jointsnvmix, we take the one that was filtered before, because checking the bam file for > 19 million variants just takes too long (~1 month)
	# /<pathTo>/filter_somatic_binom_test_JointSNVMix2 $mut_list ${outFile}
	var_dir=$out_dir/alignments_vsn_k20/variants_default/
	jointSNV=$var_dir/binomTest_160527/TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_cpp
	perc=50
	cont=20
	true_tag=""
	fn_bam=$var_dir/../bam/TU.wCont${cont}.final.RG.${perc}perc.bam
	fn_bam_normal=$var_dir/../bam/NO_final.RG.${perc}perc.bam 
	tool=$jointSNV		
	vcf_name=`basename $tool`
	vcf=$tool
	cnt_jobs_total=$((cnt_jobs_total+1))
		
	eval_dir=$var_dir/eval_${noY_}${max_cnt_normal}_binomTest_160527
	fn_out=$eval_dir/logs/`basename $vcf`_${RANDOM}.o
	mkdir -p $eval_dir
	
	if true; then
		fn_SN=$eval_dir/`basename $vcf`_indel0.eval_${eval_file_ending}
	
		exec_creation_time=`stat -c %Y $evaluate_vcf_tool_path/evaluate_vcf `
		fn_SN_creation_time=`stat -c %Y $fn_SN`	
		if [ -f $fn_SN ]; then 
			fn_SN_creation_time=`stat -c %Y $fn_SN`	
			if [ $fn_SN_creation_time -lt $exec_creation_time ]; then
				echo $fn_SN exists, and is older
			else
				echo $fn_SN exists, and is newer
				continue
			fi
		else
			echo $fn_SN does not exist
		fi
	fi
	
	command="$evaluate_vcf_tool_path/evaluate_vcf --out-dir $out_dir --vcf-dir $vcf_dir --eval-dir $eval_dir --max-cnt-normal $max_cnt_normal $tool --bam $fn_bam --bam-normal $fn_bam_normal --only-sn-auprc --pathToBam $out_dir --pathToBed $pathToBed"
	echo "$command"
	submit "$command" --tag eval_SNV --log-dir $eval_dir/logs --fn-out $fn_out --queue $queue
	cnt_jobs_TODO=$((cnt_jobs_TODO+1))
	wait_for_jobs eval_SNV --max-num $max_num_jobs
fi	


if $repeatS; then
	## Now the repeated subsampling ones
	for var_dir in $out_dir/alignments_vsn_k20/bam/variants_repeatedSubsampling/; do
		perc=50
		cont=20
		for vcf in `ls $var_dir/*.vcf | grep -v combined | grep -v paired | grep -v GERMLINE` `ls $var_dir/*.vcf.gz | grep -v combined | grep -v paired | grep -v GERMLINE`; do
			randomNum=`basename $vcf | tr '.' ' ' | awk '{print $6}' | tr '_' ' ' | awk '{print $1}'`
			fn_bam=$var_dir/../repeatedSubsampling/TU.wCont${cont}.final.RG.${perc}perc.${randomNum}.bam
			fn_bam_normal=$var_dir/../repeatedSubsampling/NO_final.RG.${perc}perc.${randomNum}.bam 
	
			cnt_jobs_total=$((cnt_jobs_total+1))
		
			eval_dir=$var_dir/eval_${noY_}${max_cnt_normal}_160527
			fn_out=$eval_dir/logs/`basename $vcf`_${RANDOM}.o			
			mkdir -p $eval_dir
	
			if true; then
				fn_SN=$eval_dir/`basename $vcf`_indel0.eval_${eval_file_ending}
	
				exec_creation_time=`stat -c %Y $evaluate_vcf_tool_path/evaluate_vcf `
				fn_SN_creation_time=`stat -c %Y $fn_SN`	
				if [ -f $fn_SN ]; then 
					fn_SN_creation_time=`stat -c %Y $fn_SN`	
					if [ $fn_SN_creation_time -lt $exec_creation_time ]; then
						echo $fn_SN exists, and is older
					else
						echo $fn_SN exists, and is newer
						continue
					fi
				else
					echo $fn_SN does not exist
				fi
			fi
	
			command="$evaluate_vcf_tool_path/evaluate_vcf --out-dir $out_dir --vcf-dir $vcf_dir --eval-dir $eval_dir --max-cnt-normal $max_cnt_normal $vcf --bam $fn_bam --bam-normal $fn_bam_normal  --only-sn-auprc --pathToBam $out_dir --pathToBed $pathToBed "
			echo "$command"
			submit "$command" --tag eval_SNV --log-dir $eval_dir/logs --fn-out $fn_out --queue $queue
			cnt_jobs_TODO=$((cnt_jobs_TODO+1))
			wait_for_jobs eval_SNV --max-num $max_num_jobs
		done
	done
	#wait_for_jobs eval_SNV --max-num $max_num_jobs
	echo "cnt_jobs_TODO=$cnt_jobs_TODO"
#exit 0
fi


#wait_for_jobs eval_SNV --max-num $max_num_jobs
echo "cnt_jobs_TODO=$cnt_jobs_TODO"
exit 0



