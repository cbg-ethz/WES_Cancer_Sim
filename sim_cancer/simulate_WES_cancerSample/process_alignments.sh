#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`
source `find ${gitDir} -name paths.sh`
source `find $SNV_DIR -name utils.sh`
source `find $SNV_DIR -name submit.sh`

bamDir=$1

## Check if sorting of bam files worked and delete otherwise
for bam in $bamDir/*.sorted.bam; do 
	if [ ! -s $bam ]; then 
		echo empty $bam
		mv $bam ${bam}_ERROR
	fi
	size_new_bam=`stat -c%s $bam`
	size_old_bam=`stat -c%s ${bam%.sorted.bam}.bam`
	size_perc=`echo "scale=3; ($size_new_bam/$size_old_bam)*100" | bc`
	if [  `echo $size_perc'<'70 | bc -l` == 1 ]; then
		echo small $bam
		mv $bam ${bam}_SMALL
	fi
done


# # sort the bam files if it failed during the alignment
for fn_bam in `find $bamDir -name *.bam | grep -v sorted | grep -v woMM | grep -v karyotypic | grep bam\/[0-9]`; do
	if [ -f ${fn_bam%.bam}.sorted.bam ]; then 
		echo "${fn_bam%.bam}.sorted.bam exists."
		continue
	fi
	command="source `find $SNV_DIR -name sort_bam.sh`; source `find $SNV_DIR -name paths.sh`; sort_bam $fn_bam $bamDir $statusDir --bam-ext .bam "
	submit "$command"  --tag sort_`basename ${fn_bam%.bam}` --log-dir $bamDir/logs/ --work-dir $SNV_DIR --queue mpi01.q
done
wait_for_jobs sort_


## Check if removing MM worked during the former run(s) of this script
for bam in $bamDir/*.woMM.bam; do 
	if [ ! -s $bam ]; then 
		echo empty $bam
		mv $bam ${bam}_ERROR
	fi
	size_new_bam=`stat -c%s $bam`
	size_old_bam=`stat -c%s ${bam%.woMM.bam}.bam`
	size_perc=`echo "scale=3; ($size_new_bam/$size_old_bam)*100" | bc` 
	echo $size_perc
	if [  `echo $size_perc'<'50 | bc -l` == 1 ]; then
		echo small $bam
		mv $bam ${bam}_SMALL
	fi
done

# echo remove multimappers
cnt=0
log_dir=$bamDir/logs/
mkdir -p $log_dir
for bam_file in $bamDir/*.sorted.bam
do
	fn_out=${bam_file%.bam}.woMM.bam
	if [ -f $fn_out ]; then
		echo "$fn_out already exists."
		continue;
	fi
	command="$samtools view -bh -F 256 $bam_file > $fn_out" 
	echo "submit $command  --tag rmMM`basename ${bam_file%.bam}` --log-dir $log_dir"
	submit "$command"  --tag rmMM`basename ${bam_file%.bam}` --log-dir $log_dir --queue mpi01.q
done
wait_for_jobs rmMM
	

## check if some files had an error while sorting karyotypically during a former run(s) of this script
for bam in `find $bamDir -name *karyotypic.bam`; do 
	if [ ! -s $bam ]; then 
		echo "empty:"
		echo $bam
		mv $bam ${bam}_ERROR
	fi
	size_new_bam=`stat -c%s $bam`
	size_old_bam=`stat -c%s ${bam%.karyotypic.bam}.bam`
	size_perc=`echo "scale=3; ($size_new_bam/$size_old_bam)*100" | bc` 
	echo $size_perc
	if [  `echo $size_perc'<'70 | bc -l` == 1 ]; then
		echo small $bam
		mv $bam ${bam}_SMALL
	fi
done

# # sort the bam files (already done after alignment in our pipeline)
# echo sort them in karyotypic order
for fn_bam in `find $bamDir -name *.sorted.woMM.bam`; do
	if [ -f ${fn_bam%.bam}.karyotypic.bam ]; then 
		continue
	fi
	command="source `find $SNV_DIR -name sort_bam.sh`; source `find $SNV_DIR -name paths.sh`; sort_bam $fn_bam $bamDir $statusDir  --bam-ext .sorted.woMM.bam --human-karyotypic"
	submit "$command"  --tag sort_`basename ${fn_bam%.bam}` --log-dir $bamDir/logs/ --work-dir $SNV_DIR --queue mpi01.q
done
wait_for_jobs sort_


# echo merge the TU and NO bam files together
outBamTU=$bamDir/TU_final.bam
if [ -f $outBamTU ]; then
	echo "$outBamTU already exists."
else
	command="$samtools merge $outBamTU $bamDir/*TU*.sorted.woMM.karyotypic.bam"
	submit "$command"  --tag merge_job_`basename ${outBamTU%.bam}` --log-dir $bamDir/logs/ --queue mpi01.q
fi

outBamNO=$bamDir/NO_final.bam
if [ -f $outBamNO ]; then
		echo "$outBamNO already exists."
else
	command="$samtools merge $outBamNO $bamDir/*NO*.sorted.woMM.karyotypic.bam"
	submit "$command"  --tag merge_job_`basename ${outBamNO%.bam}` --log-dir $bamDir/logs/ --queue mpi01.q
fi
wait_for_jobs merge_job_


## echo replace read groups such that we only have one read group in the header (with only one sample name)
for curr_bam in $bamDir/NO_final.bam
do
	echo $curr_bam
	outputBam=${curr_bam%.bam}.RG.bam

	if [ -f $outputBam ]; then 
		continue
	fi

	RGID=`basename ${curr_bam%.bam}`
	RGLB=1
	RGPL="illumina"
	RGPU=1
	RGSM=`basename ${curr_bam%.bam} | tr '.' ' '  |  tr '_' ' ' | awk '{print $1}'`

	command="source `find $SNV_DIR -name replaceReadGroups.sh`; source `find $SNV_DIR -name paths.sh`; replaceReadGroups $curr_bam $outputBam $RGID $RGLB $RGPL $RGPU $RGSM --check-fexist"
	submit "$command"  --tag replace_RG`basename ${curr_bam%.bam}` --log-dir $bamDir/logs/  --queue mpi01.q
done
wait_for_jobs replace_RG


echo create subsets of the bam files to merge with the TU file
bamFile=$bamDir/NO_final.bam
for fraction in 0.10 0.20 0.40 0.60; do 
	outBam=${bamFile%.bam}.`echo $fraction | cut -f 2 -d "."`perc.bam
	if [ -f $outBam ]; then 
		continue; 
	fi
	command="source `find $SNV_DIR -name subsetSam.sh`; source `find $SNV_DIR -name paths.sh`; subsetSam $bamFile $outBam $fraction"
	submit "$command"  --tag subset1NO --log-dir $bamDir/logs/ --queue mpi01.q
done
wait_for_jobs subset1

# echo merge the TU final bam in the bamDir directory with the 10%, 20%, 40%, and 60% bam file
# to simulate normal contamination
for cont in 10 20 40 60
do
	outBamMerge=$bamDir/TU.wCont${cont}.final.bam
	if [ -f $outBamMerge ]; then
		echo "$outBamMerge already exists."
	else
		echo $bamFilesMerge
		command="$samtools merge $outBamMerge $bamDir/TU_final.bam $bamDir/NO_final.${cont}perc.bam"
		submit "$command"  --tag merge_TU_`basename ${outBamMerge%.bam}` --log-dir $bamDir/logs/ --queue mpi01.q
	fi
done
wait_for_jobs merge_TU_

## echo replace read groups such that we only have one read group in the header (with only one sample name)
for curr_bam in $bamDir/TU.wCont[0-9]0.final.bam
do
	echo $curr_bam
	outputBam=${curr_bam%.bam}.RG.bam

	if [ -f $outputBam ]; then 
		continue
	fi

	RGID=`basename ${curr_bam%.bam}`
	RGLB=1
	RGPL="illumina"
	RGPU=1
	RGSM=`basename ${curr_bam%.bam} | tr '.' ' '  |  tr '_' ' ' | awk '{print $1}'`

	command="source `find $SNV_DIR -name replaceReadGroups.sh`; source `find $SNV_DIR -name paths.sh`; replaceReadGroups $curr_bam $outputBam $RGID $RGLB $RGPL $RGPU $RGSM --check-fexist"
	submit "$command"  --tag replace_RG`basename ${curr_bam%.bam}` --log-dir $bamDir/logs/  --queue mpi01.q
done
wait_for_jobs replace_RG

      

for frac in 12 25 50 75 100; do

	all_cont=20
	if [ $frac -eq 50 ]; then 
		all_cont="10 20 40 60"
		echo run for all: $all_cont
	fi
	for cont in $all_cont; do
		bamFile=$bamDir/TU.wCont${cont}.final.RG.bam
		fraction=`echo "scale=6; $frac/1.${cont}/100" | bc`
		outBam=${bamFile%.bam}.${frac}perc.bam
		echo create file $outBam
		if [ -f $outBam ]; then 
			continue; 
		fi
		command="source `find $SNV_DIR -name subsetSam.sh`; source `find $SNV_DIR -name paths.sh`; subsetSam $bamFile $outBam $fraction; $samtools index $outBam"
		submit "$command"  --tag subsetTU --log-dir $bamDir/logs/ --queue mpi01.q
	done

	#Normal files
	bamFile=$bamDir/NO_final.RG.bam
	fraction=`echo "scale=6; $frac/100" | bc`
	outBam=${bamFile%.bam}.${frac}perc.bam
	if [ -f $outBam ]; then 
		continue; 
	fi
	command="source `find $SNV_DIR -name subsetSam.sh`; source `find $SNV_DIR -name paths.sh`; subsetSam $bamFile $outBam $fraction; $samtools index $outBam"
	submit "$command"  --tag subsetNO --log-dir $bamDir/logs/ --queue mpi01.q

done
wait_for_jobs subset

