#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`

varDir=$base_dir/alignments_vsn_k20/variants_default/
evalDir=$varDir/eval_10000_160527

outDir=$varDir/all_pairwise_overlaps/
mkdir -p $outDir
binomDir=$varDir/eval_10000_binomTest_160527
tuneDir=$(echo $varDir | sed 's/\_default/\_tuned/')
#echo $tuneDir
for vcf in $varDir/TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf $tuneDir/TU.wCont20.final.RG.50perc.bam_minvarfreq0.02_varscan2.txt.snp.Somatic_qual.vcf $binomDir/TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_binomTestSomatic.vcf $tuneDir/TU.wCont20.final.RG.50perc.option_C200_noE_samvar_1_2.SOMATIC.vcf $varDir/TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf $binomDir/TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_binomTestSomatic.vcf $varDir/TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf ${tuneDir}/sinvictqscorecutoff60_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf $varDir/TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf; do
	#ls -lh $vcf
	if [ ! -f $outDir/`basename $vcf` ]; then
		ln -s $vcf $outDir/`basename $vcf`
	fi
done

# overlaps
mkdir -p $outDir/pairwise_overlaps

cd $outDir/pairwise_overlaps
allVCFs=""
cnt=0
for vcf in ../*.vcf; do
	echo $vcf
	sorted_pos=$(echo `basename ${vcf%.vcf} | sed 's/TU.wCont20.final.RG.50perc.//' | sed 's/TU.wCont20.final.RG.50perc_//' | sed 's/SOMATIC//' | sed 's/noComma//' | sed 's/Raw//' | sed 's/rewritten//' | sed 's/txt.snp.Somatic_qual//' | sed 's/binomTestSomatic//' | sed 's/NO_final.RG.50perc_alternative_two.sided.//' | sed 's/option_C200_noE_//' | sed 's/SNVs.raw..//' | sed 's/vs_NO_final.RG.50perc_//' `_pos.sorted.txt)
	if [ ! -f $sorted_pos ]; then
		cat $vcf | grep -v ^# | sort -k1,1 -k2n | awk 'BEGIN { OFS = "\t"; ORS = "\n" } {print $1,$2}' | uniq > $sorted_pos
	else
		echo "$sorted_pos already exists"
	fi
	allVCFs[$cnt]=$sorted_pos
	cnt=$((cnt+1))
done

echo "These were all files found for pairwise overlap: $allVCFs[*]"

module load java/jdk7
for i in `seq 1 9`; do
	echo ${allVCFs[i-1]}
done


# all pairwise overlaps - there are 36
for numIdx in 01 02 03 04 05 06 07 08 12 13 14 15 16 17 18 23 24 25 26 27 28 34 35 36 37 38 45 46 47 48 56 57 58 67 68 78; do
	idxIn1=`echo $numIdx | grep -o . | awk 'NR==1'`
        idxIn2=`echo $numIdx | grep -o . | awk 'NR==2'`
		
	cnt=0
	file1=${allVCFs[$idxIn1]}
	file2=${allVCFs[$idxIn2]}
	overlapFile=overlap_`basename ${file1%.txt}`_`basename ${file2}`
	privateFile=`echo $overlapFile | sed 's/overlap/private/'`
	if [ ! -f $overlapFile ]; then
		command="java SNV_Overlapper_Ordered $file1 $file2 $overlapFile $privateFile"
		echo $commmand
		eval $command
	fi
	cat $overlapFile | tr '\t' '_' > ${overlapFile%.txt}_.txt
	if [[ $file1 == *joint* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf
	elif [[ $file1 == *deepSNV* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_binomTestSomatic.vcf
	elif [[ $file1 == *varscan2* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.bam_minvarfreq0.02_varscan2.txt.snp.Somatic_qual.vcf
	elif [[ $file1 == *muTect* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_binomTestSomatic.vcf
	elif [[ $file1 == *samvar* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.option_C200_noE_samvar_1_2.SOMATIC.vcf
	elif [[ $file1 == *gatk_SNVS* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf
	elif [[ $file1 == *gatkHPCaller* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf
	elif [[ $file1 == *somaticSniper* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf
	elif [[ $file1 == *sinvict* ]]; then
		originalVCF=../sinvictqscorecutoff60_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf
	fi
	overlapFileVCForiginal=${overlapFile%.txt}_original.vcf
	if [ ! -f $overlapFileVCForiginal ]; then
		command="java createOriginalVCF ${overlapFile%.txt}_.txt $originalVCF $overlapFileVCForiginal"
		echo $command
		eval $command
	fi

	if [[ $file2 == *joint* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf
	elif [[ $file2 == *deepSNV* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_binomTestSomatic.vcf
	elif [[ $file2 == *varscan2* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.bam_minvarfreq0.02_varscan2.txt.snp.Somatic_qual.vcf
	elif [[ $file2 == *muTect* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_binomTestSomatic.vcf
	elif [[ $file2 == *samvar* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.option_C200_noE_samvar_1_2.SOMATIC.vcf
	elif [[ $file2 == *gatk_SNVS* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf
	elif [[ $file2 == *gatkHPCaller* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf
	elif [[ $file2 == *somaticSniper* ]]; then
		originalVCF=../TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf
	elif [[ $file2 == *sinvict* ]]; then
		originalVCF=../sinvictqscorecutoff60_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf
	fi
	overlapFileVCForiginal=${overlapFile%.txt}_original_file2.vcf
	if [ ! -f $overlapFileVCForiginal ]; then
		command="java createOriginalVCF ${overlapFile%.txt}_.txt $originalVCF $overlapFileVCForiginal"
		echo $command
		eval $command
	fi
done
## print all subsets of size two of all nine callers:
declare -a allTools=('deepSNV' 'GATK HP' 'GATK UG' 'JointSNVMix2' 'MuTect' 'SAMtools' 'SiNVICT' 'somaticSniper' 'VarScan2')
declare -a allToolsnamesinfiles=('deepSNV' 'gatkHP' 'gatk_SNVs' 'joint' 'muTect' 'samvar' 'sinvict' 'somaticSniper' 'varscan2' )
for numIdx in 01 02 03 04 05 06 07 08 12 13 14 15 16 17 18 23 24 25 26 27 28 34 35 36 37 38 45 46 47 48 56 57 58 67 68 78; do 
	idxIn1=`echo $numIdx | grep -o . | awk 'NR==1'`; 
	idxIn2=`echo $numIdx | grep -o . | awk 'NR==2'`; 
	tool1=${allToolsnamesinfiles[$idxIn1]}
	tool2=${allToolsnamesinfiles[$idxIn2]}
	overlapfile=`ls overlap*_pos.sorted.txt | grep $tool1 | grep $tool2`	
	#echo "overlapfile = $overlapfile"
	numSharedCalls=`cat $overlapfile | wc -l`
	file1forRankingauPRC10File=`ls ./eval_10000_160527/overlap*original.vcf_indel0.eval_auPRC_10 | grep $tool1 | grep $tool2`
	#echo "file1forRankingauPRC10File = $file1forRankingauPRC10File"
	file2forRankingauPRC10File=`ls ./eval_10000_160527/overlap*original_file2.vcf_indel0.eval_auPRC_10 | grep $tool1 | grep $tool2`
	#echo "file2forRankingauPRC10File = $file2forRankingauPRC10File"
	file1rankingauPRC=`cat $file1forRankingauPRC10File`
	if [ ! -f $file2forRankingauPRC10File ]; then
		file2rankingauPRC=0
	else
		file2rankingauPRC=`cat $file2forRankingauPRC10File`
	fi
	fileused1=`basename $file1forRankingauPRC10File | sed 's/_indel0.eval_auPRC_10//' | sed 's/overlap_//'`
	fileused2=`basename $file2forRankingauPRC10File | sed 's/_indel0.eval_auPRC_10//' | sed 's/overlap_//'`
	if [ `echo $file1rankingauPRC'>='$file2rankingauPRC | bc -l` == 1 ]; then 
		echo "${allTools[$idxIn1]}, ${allTools[$idxIn2]} & $numSharedCalls & $fileused1  & \num{$file1rankingauPRC} "; 
	else
		echo "${allTools[$idxIn1]}, ${allTools[$idxIn2]} & $numSharedCalls & $fileused2  &  \num{$file2rankingauPRC}  ";
	fi
done

cd -

