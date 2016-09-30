#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`

varDir=$base_dir/alignments_vsn_k20/variants_default/
evalDir=$varDir/eval_10000_160527
specific_overlaps=false
perc5=false
perc1=false
perc10=true


if $specific_overlaps; then
	## only specific overlaps:
	outDir=$varDir/specific_overlaps
	mkdir -p $outDir
	binomDir=$varDir/eval_10000_binomTest_160527
	tuneDir=$(echo $varDir | sed 's/\_default/\_tuned/')
	#echo $tuneDir
	for vcf in $varDir/TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf $binomDir/TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_binomTestSomatic.vcf $tuneDir/TU.wCont20.final.RG.50perc.bam_minvarfreq0.02_varscan2.txt.snp.Somatic_qual.vcf $binomDir/TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_binomTestSomatic.vcf $tuneDir/TU.wCont20.final.RG.50perc.option_C200_noE_samvar_1_2.SOMATIC.vcf $varDir/TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf $varDir/TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf $varDir/TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf ${tuneDir}/sinvictqscorecutoff60_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf; do
		#ls -lh $vcf
		if [ ! -f $outDir/`basename $vcf` ]; then
			ln -s $vcf $outDir/`basename $vcf`
		fi
	done

	# overlaps
	mkdir -p $outDir/overlaps

	cd $outDir/overlaps
		for vcf in ../*.vcf; do
			echo $vcf
			sorted_pos=$(echo `basename ${vcf%.vcf} | sed 's/TU.wCont20.final.RG.50perc.//' | sed 's/TU.wCont20.final.RG.50perc_//' | sed 's/SOMATIC//' | sed 's/noComma//' | sed 's/Raw//' | sed 's/rewritten//' | sed 's/txt.snp.Somatic_qual//' | sed 's/binomTestSomatic//' | sed 's/NO_final.RG.50perc_alternative_two.sided.//' | sed 's/option_C200_noE_//' | sed 's/SNVs.raw..//' | sed 's/vs_NO_final.RG.50perc_//' `_pos.sorted.txt)
			if [ ! -f $sorted_pos ]; then
				cat $vcf | grep -v ^# | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > $sorted_pos
			else
				echo "$sorted_pos already exists"
			fi
		done

	module load java/jdk7
	file1=deepSNV.vcf__pos.sorted.txt
	file2=jointSNVMix2_SNVs__pos.sorted.txt
	overlapFile=overlap_`basename ${file1%.txt}`_`basename ${file2}`
	privateFile=private_`basename ${file1%.txt}`_`basename ${file2}`
	if [ ! -f $overlapFile ]; then
		command="java SNV_Overlapper_Ordered $file1 $file2 $overlapFile $privateFile"
		echo $command
		eval $command
	fi

	file1=overlap_deepSNV.vcf__pos.sorted_jointSNVMix2_SNVs__pos.sorted.txt
	file2=bam_minvarfreq0.02_varscan2._pos.sorted.txt
	overlapFile=overlap_`basename ${file1%.txt}`_`basename ${file2}`
	privateFile=private_`basename ${file1%.txt}`_`basename ${file2}`
	if [ ! -f $overlapFile ]; then
		command="java SNV_Overlapper_Ordered $file1 $file2 $overlapFile $privateFile"
		echo $command
		eval $command
	fi

	file1=overlap_overlap_deepSNV.vcf__pos.sorted_jointSNVMix2_SNVs__pos.sorted_bam_minvarfreq0.02_varscan2._pos.sorted.txt
	file2=muTect_SNVs_.vcf__pos.sorted.txt
	overlapFile=overlap_deepSNV_joint_varscan_muTect.txt
	privateFile=private_deepSNV_joint_varscan_muTect.txt
	if [ ! -f $overlapFile ]; then
		command="java SNV_Overlapper_Ordered $file1 $file2 $overlapFile $privateFile"
		echo $command
		eval $command
	fi
	#cat $overlapFile | tr '\t' '_' > ${overlapFile%.txt}_.txt
	#command="java createOriginalVCF ${overlapFile%.txt}_.txt ../TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf ${overlapFile%.txt}_original_joint.vcf"
	#echo $command
	#eval $command

	file1=overlap_deepSNV_joint_varscan_muTect.txt
	file2=samvar_1_2._pos.sorted.txt
	overlapFile=overlap_deepSNV_joint_varscan_muTect_samvar.txt
	privateFile=`echo $overlapFile | sed 's/overlap/private/'`
	if [ ! -f $overlapFile ]; then
		command="java SNV_Overlapper_Ordered $file1 $file2 $overlapFile $privateFile"
		echo $command
		eval $command
	fi
	cat $overlapFile | tr '\t' '_' > ${overlapFile%.txt}_.txt
	command="java createOriginalVCF ${overlapFile%.txt}_.txt ../TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf ${overlapFile%.txt}_original_joint.vcf; mv ${overlapFile%.txt}_original_joint.vcf `echo ${overlapFile%.txt}_original_joint.vcf | sed 's/samvar/samvar_1_2/'`"
	echo $command
	eval $command

	file1=overlap_deepSNV_joint_varscan_muTect.txt
	file2=somaticSniper_SNVs__qual._pos.sorted.txt
	#overlapFile=overlap_`basename ${file1%.txt}`_`basename ${file2}`
	#privateFile=private_`basename ${file1%.txt}`_`basename ${file2}`
	overlapFile=overlap_deepSNV_joint_varscan_muTect_somaticSniper.txt
	privateFile=`echo $overlapFile | sed 's/overlap/private/'`
	if [ ! -f $overlapFile ]; then
		command="java SNV_Overlapper_Ordered $file1 $file2 $overlapFile $privateFile"
		echo $command
		eval $command
	fi
	cat $overlapFile | tr '\t' '_' > ${overlapFile%.txt}_.txt
	command="java createOriginalVCF ${overlapFile%.txt}_.txt ../TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf ${overlapFile%.txt}_original_joint.vcf"
	echo $command
	eval $command

	file1=overlap_deepSNV_joint_varscan_muTect.txt
	file2=sinvictqscorecutoff60_somatic_pos.sorted.txt
	overlapFile=overlap_deepSNV_joint_varscan_muTect_sinvict.txt
	privateFile=`echo $overlapFile | sed 's/overlap/private/'`
	if [ ! -f $overlapFile ]; then
		command="java SNV_Overlapper_Ordered $file1 $file2 $overlapFile $privateFile"
		echo $command
		eval $command
	fi
	cat $overlapFile | tr '\t' '_' > ${overlapFile%.txt}_.txt
	command="java createOriginalVCF ${overlapFile%.txt}_.txt ../TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf ${overlapFile%.txt}_original_joint.vcf"
	echo $command
	eval $command


	file1=overlap_deepSNV_joint_varscan_muTect_sinvict.txt
	file2=somaticSniper_SNVs__qual._pos.sorted.txt
	overlapFile=overlap_`basename ${file1%.txt}`_`basename ${file2}`
	privateFile=private_`basename ${file1%.txt}`_`basename ${file2}`
	if [ ! -f $overlapFile ]; then
		command="java SNV_Overlapper_Ordered $file1 $file2 $overlapFile $privateFile"
		echo $command
		eval $command
	fi

	file1=overlap_overlap_deepSNV_joint_varscan_muTect_sinvict_somaticSniper_SNVs__qual._pos.sorted.txt
	file2=gatk_SNVs.raw._pos.sorted.txt
	overlapFile=overlap_`basename ${file1%.txt}`_`basename ${file2}`
	privateFile=private_`basename ${file1%.txt}`_`basename ${file2}`
	if [ ! -f $overlapFile ]; then
		command="java SNV_Overlapper_Ordered $file1 $file2 $overlapFile $privateFile"
		echo $command
		eval $command
	fi

	file1=overlap_overlap_overlap_deepSNV_joint_varscan_muTect_sinvict_somaticSniper_SNVs__qual._pos.sorted_gatk_SNVs.raw._pos.sorted.txt
	file2=gatkHPCaller__pos.sorted.txt
	#overlapFile=overlap_`basename ${file1%.txt}`_`basename ${file2}`
	#privateFile=private_`basename ${file1%.txt}`_`basename ${file2}`
	overlapFile=overlap_deepSNV_joint_varscan_muTect_sinvict_somaticSniper_gatk_SNVs_gatkHPCaller.txt
	privateFile=`echo $overlapFile | sed 's/overlap/private/'`
	if [ ! -f $overlapFile ]; then
		command="java SNV_Overlapper_Ordered $file1 $file2 $overlapFile $privateFile"
		echo $command
		eval $command
	fi

	file1=overlap_deepSNV_joint_varscan_muTect_sinvict_somaticSniper_gatk_SNVs_gatkHPCaller.txt
	file2=samvar_1_2._pos.sorted.txt
	#overlapFile=overlap_`basename ${file1%.txt}`_`basename ${file2}`
	#privateFile=private_`basename ${file1%.txt}`_`basename ${file2}`
	overlapFile=overlap_deepSNV_joint_varscan_muTect_sinvict_somaticSniper_gatk_SNVs_gatkHPCaller_samvar.txt
	privateFile=`echo $overlapFile | sed 's/overlap/private/'`
	if [ ! -f $overlapFile ]; then
		command="java SNV_Overlapper_Ordered $file1 $file2 $overlapFile $privateFile"
		echo $command
		eval $command
	fi
	cat $overlapFile | tr '\t' '_' > ${overlapFile%.txt}_.txt
	command="java createOriginalVCF ${overlapFile%.txt}_.txt ../TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf ${overlapFile%.txt}_original_joint.vcf; mv ${overlapFile%.txt}_original_joint.vcf `echo ${overlapFile%.txt}_original_joint.vcf | sed 's/samvar/samvar_1_2/'`"
	echo $command
	eval $command

	cd -
exit 0
fi


if $perc5; then
	cd $varDir
	mkdir -p SNV_lists_fdr05_5percent
	
	for snv_list in $evalDir/TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_fdr05.vcf $evalDir/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_fdr05.vcf $evalDir/TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_fdr05.vcf $evalDir/TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_fdr05.vcf $evalDir/TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_fdr05.vcf
	do
		if [ ! -f ./SNV_lists_fdr05_5percent/`basename $snv_list` ]; then
			ln -s $snv_list ./SNV_lists_fdr05_5percent
		fi
	done
	
	
	# overlaps
	mkdir -p ./SNV_lists_fdr05_5percent/overlaps
	cd ./SNV_lists_fdr05_5percent/overlaps
	if [ -f ../TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_fdr05.vcf -a ! -f varscan2.txt ]; then
		cat ../TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_fdr05.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}'  | uniq > varscan2.txt
	fi
	if [ -f ../sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_fdr05.vcf -a ! -f sinvict.txt ]; then
		cat ../sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_fdr05.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}'  | uniq > sinvict.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_fdr05.vcf -a ! -f gatkUG.txt ]; then
		cat ../TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_fdr05.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > gatkUG.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_fdr05.vcf -a ! -f deep.txt ]; then
		# have to make them unique, because deepSNV reports quite some double SNVs at the same position -> to not screw up the overlap counts, I make them unique
		cat ../TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_fdr05.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > deep.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_fdr05.vcf -a ! -f joint.txt ]; then
		cat ../TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_fdr05.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > joint.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_fdr05.vcf -a ! -f somSniper.txt ]; then
		cat ../TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_fdr05.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > somSniper.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_fdr05.vcf -a ! -f muTect.txt ]; then
		cat ../TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_fdr05.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > muTect.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_fdr05.vcf -a ! -f gatkHP.txt ]; then
		cat ../TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_fdr05.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > gatkHP.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_fdr05.vcf -a ! -f SAMvar.txt ]; then
		cat ../TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_fdr05.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > SAMvar.txt
	fi
	cd -
	
	module load java/jdk7
	cnt=0
	for file1 in ./SNV_lists_fdr05_5percent/overlaps/*.txt; do 
		for file2 in ./SNV_lists_fdr05_5percent/overlaps/*.txt; do
			if [ ! -f ./SNV_lists_fdr05_5percent/overlaps/overlap_`basename ${file1%.txt}`_`basename ${file2}` ]; then
				if [[ $file1 != *private* && $file2 != *private*  ]]; then
					command="java SNV_Overlapper_Ordered $file1 $file2 ./SNV_lists_fdr05_5percent/overlaps/overlap_`basename ${file1%.txt}`_`basename ${file2}` ./SNV_lists_fdr05_5percent/overlaps/private_`basename ${file1%.txt}`_`basename ${file2}`"
					echo $command
					eval $command
				cnt=$((cnt+1))
				fi
			fi
			if [ $cnt -gt 50 ]; then
				break
			fi
		done
	done

	# summary file with overlaps for the top 5 callers
	summary_overlaps=./SNV_lists_fdr05_5percent/overlaps/overlap_summary_top5callers.txt
	if [ ! -f $summary_overlaps ]; then
		${gitDir}/sim_cancer/combis_overlaps/extract_nums_forVennDiagram.sh ./SNV_lists_fdr05_5percent/overlaps/ > $summary_overlaps
		${gitDir}/sim_cancer/combis_overlaps/extract_nums_forVennDiagram.sh ./SNV_lists_fdr05_5percent/overlaps/ > $summary_overlaps
	fi
	
	cd $currScriptDir
exit 0
fi




if $perc1; then
	cd $varDir
	mkdir -p SNV_lists_fdr01_1percent
	
	for snv_list in  $evalDir/TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_fdr01.vcf $evalDir/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_fdr01.vcf  $evalDir/TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_fdr01.vcf $evalDir/TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_fdr01.vcf $evalDir/TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_fdr01.vcf
	do
		if [ ! -f ./SNV_lists_fdr01_1percent/`basename $snv_list` ]; then
			ln -s $snv_list ./SNV_lists_fdr01_1percent
		fi
	done
	
	# overlaps
	mkdir -p ./SNV_lists_fdr01_1percent/overlaps
	cd ./SNV_lists_fdr01_1percent/overlaps
	if [ -f ../TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_fdr01.vcf -a ! -f varscan2.txt ]; then
		cat ../TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_fdr01.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > varscan2.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_fdr01.vcf -a ! -f gatkUG.txt ]; then
		#cat ../TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC_paired.noComma.vcf_indel0.eval_SN_fdr01.vcf | sort -k1,1 -k2n > gatkUG.txt
		cat ../TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_fdr01.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > gatkUG.txt
	fi
	if [ -f ../sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_fdr01.vcf -a ! -f sinvict.txt ]; then
		cat ../sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_fdr01.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > sinvict.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_fdr01.vcf -a ! -f deep.txt ]; then
		# have to make them unique, because deepSNV reports quite some double SNVs at the same position -> to not screw up the overlap counts, I make them unique
		cat ../TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_fdr01.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > deep.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_fdr01.vcf -a ! -f joint.txt ]; then
		cat ../TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_fdr01.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > joint.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_fdr01.vcf -a ! -f somSniper.txt ]; then
		#cat ../TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_fdr01.vcf | sort -k1,1 -k2n > somSniper.txt
		cat ../TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_fdr01.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > somSniper.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_fdr01.vcf -a ! -f muTect.txt ]; then
		cat ../TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_fdr01.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > muTect.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_fdr01.vcf -a ! -f gatkHP.txt ]; then
		#cat ../TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_fdr01.vcf | sort -k1,1 -k2n > gatkHP.txt
		cat ../TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_fdr01.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > gatkHP.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_fdr01.vcf -a ! -f SAMvar.txt ]; then
		#cat ../TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_fdr01.vcf | sort -k1,1 -k2n > SAMvar.txt
		cat ../TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_fdr01.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > SAMvar.txt
	fi
	cd -
	
	module load java/jdk7
	cnt=0
	for file1 in ./SNV_lists_fdr01_1percent/overlaps/*.txt; do 
		for file2 in ./SNV_lists_fdr01_1percent/overlaps/*.txt; do
			if [ ! -f ./SNV_lists_fdr01_1percent/overlaps/overlap_`basename ${file1%.txt}`_`basename ${file2}` ]; then
				if [[ $file1 != *private* && $file2 != *private*  ]]; then
					command="java SNV_Overlapper_Ordered $file1 $file2 ./SNV_lists_fdr01_1percent/overlaps/overlap_`basename ${file1%.txt}`_`basename ${file2}` ./SNV_lists_fdr01_1percent/overlaps/private_`basename ${file1%.txt}`_`basename ${file2}`"
					echo $command
					eval $command
				fi
			fi
			cnt=$((cnt+1))
			if [ $cnt -gt 100 ]; then
				break;
			fi
		done
		if [ $cnt -gt 100 ]; then
			break
		fi
	done


	# summary file with overlaps for the top 5 callers
	summary_overlaps=./SNV_lists_fdr01_1percent/overlaps/overlap_summary_top5callers.txt
	if [ ! -f $summary_overlaps ]; then
		${gitDir}/sim_cancer/combis_overlaps/extract_nums_forVennDiagram.sh ./SNV_lists_fdr01_1percent/overlaps/ > $summary_overlaps
		${gitDir}/sim_cancer/combis_overlaps/extract_nums_forVennDiagram.sh ./SNV_lists_fdr01_1percent/overlaps/ > $summary_overlaps
	fi
	cd $currScriptDir
exit 0
fi


if $perc10; then
	cd $varDir
	mkdir -p SNV_lists_fdr10_10percent
	
	for snv_list in $evalDir/TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_fdr10.vcf  $evalDir/sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_fdr10.vcf  $evalDir/TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_fdr10.vcf $evalDir/TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_fdr10.vcf $evalDir/TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_fdr10.vcf
	do
		if [ ! -f ./SNV_lists_fdr10_10percent/`basename $snv_list` ]; then
			ln -s $snv_list ./SNV_lists_fdr10_10percent
		fi
	done
	
	# overlaps
	mkdir -p ./SNV_lists_fdr10_10percent/overlaps
	cd ./SNV_lists_fdr10_10percent/overlaps
	if [ -f ../TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_fdr10.vcf -a ! -f varscan2.txt ]; then
		cat ../TU.wCont20.final.RG.50perc.bam__varscan2.txt.snp.Somatic_qual.vcf_indel0.eval_SN_fdr10.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > varscan2.txt
	fi
	if [ -f ../sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_fdr10.vcf -a ! -f sinvict.txt ]; then
		cat ../sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf_indel0.eval_SN_fdr10.vcf  | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > sinvict.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_fdr10.vcf -a ! -f gatkUG.txt ]; then
		#cat ../TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC_paired.noComma.vcf_indel0.eval_SN_fdr10.vcf | sort -k1,1 -k2n > gatkUG.txt
		cat ../TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf_indel0.eval_SN_fdr10.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > gatkUG.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_fdr10.vcf -a ! -f deep.txt ]; then
		# have to make them unique, because deepSNV reports quite some double SNVs at the same position -> to not screw up the overlap counts, I make them unique
		cat ../TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_indel0.eval_SN_fdr10.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > deep.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_fdr10.vcf -a ! -f joint.txt ]; then
		cat ../TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf_indel0.eval_SN_fdr10.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > joint.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_fdr10.vcf -a ! -f somSniper.txt ]; then
		#cat ../TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_fdr10.vcf | sort -k1,1 -k2n > somSniper.txt
		cat ../TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf_indel0.eval_SN_fdr10.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > somSniper.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_fdr10.vcf -a ! -f muTect.txt ]; then
		cat ../TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_indel0.eval_SN_fdr10.vcf | sort -k1,1 -k2n > muTect.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_fdr10.vcf -a ! -f gatkHP.txt ]; then
		#cat ../TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_fdr10.vcf | sort -k1,1 -k2n > gatkHP.txt
		cat ../TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf_indel0.eval_SN_fdr10.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > gatkHP.txt
	fi
	if [ -f ../TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_fdr10.vcf -a ! -f SAMvar.txt ]; then
		#cat ../TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_fdr10.vcf | sort -k1,1 -k2n > SAMvar.txt
		cat ../TU.wCont20.final.RG.50perc.option__noE_samvar_1_2.SOMATIC.vcf.gz_indel0.eval_SN_fdr10.vcf | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > SAMvar.txt
	fi
	cd -
	
	module load java/jdk7
	cnt=0
	for file1 in ./SNV_lists_fdr10_10percent/overlaps/*.txt; do 
		for file2 in ./SNV_lists_fdr10_10percent/overlaps/*.txt; do
			if [ ! -f ./SNV_lists_fdr10_10percent/overlaps/overlap_`basename ${file1%.txt}`_`basename ${file2}` ]; then
				if [[ $file1 != *private* && $file2 != *private*  ]]; then
					command="java SNV_Overlapper_Ordered $file1 $file2 ./SNV_lists_fdr10_10percent/overlaps/overlap_`basename ${file1%.txt}`_`basename ${file2}` ./SNV_lists_fdr10_10percent/overlaps/private_`basename ${file1%.txt}`_`basename ${file2}`"
					echo $command
					eval $command
				fi
			fi
			cnt=$((cnt+1))
			if [ $cnt -gt 100 ]; then
				break;
			fi
		done
		if [ $cnt -gt 100 ]; then
			break
		fi
	done

	# summary file with overlaps for the top 5 callers
	summary_overlaps=./SNV_lists_fdr10_10percent/overlaps/overlap_summary_top5callers.txt
	if [ ! -f $summary_overlaps ]; then
		${gitDir}/sim_cancer/combis_overlaps/extract_nums_forVennDiagram.sh ./SNV_lists_fdr10_10percent/overlaps/ > $summary_overlaps
		${gitDir}/sim_cancer/combis_overlaps/extract_nums_forVennDiagram.sh ./SNV_lists_fdr10_10percent/overlaps/ > $summary_overlaps
	fi
	cd $currScriptDir
exit 0
fi
