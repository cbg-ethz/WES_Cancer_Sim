#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`
source `find ${gitDir} -name paths.sh`
source `find $SNV_DIR -name utils.sh`
source `find $SNV_DIR -name submit.sh`

bamDir=$1
bedDir=$2

if [ -z "$1" -o -z "$2" ]; then
	echo "Please provide bamDir and bedDir"
	exit 0
fi


## Check if IDsorting of bam files worked and delete otherwise
for bam in $bamDir/*.IDsorted.bam; do 
	if [ ! -f $bam ]; then
		echo "No IDsorted bam file exists yet"	
		continue
	fi
	if [ ! -s $bam ]; then 
		echo empty $bam
		mv $bam ${bam}_ERROR
	fi
	size_new_bam=`stat -c%s $bam`
	size_old_bam=`stat -c%s ${bam%.IDsorted.bam}.bam`
	size_perc=`echo "scale=3; ($size_new_bam/$size_old_bam)*100" | bc`
	if [  `echo $size_perc'<'70 | bc -l` == 1 ]; then
		echo small $bam
		mv $bam ${bam}_SMALL
	fi
done


for f in $bamDir/*.sorted.bam; do
	fn_IDsorted=${f%.bam}.IDsorted.bam
	if [ ! -f $fn_IDsorted ]; then 
		fn_out=$bamDir/logs/`basename ${f%.bam}`_IDsort.o
		command="${samtools_1_2} sort -n -l 9 $f ${fn_IDsorted%.bam}"
		submit "$command" --tag IDsort --work-dir `pwd` --log-dir $bamDir/logs/ --queue mpi01.q --fn-out $fn_out
	fi
done
wait_for_jobs IDsort

for f in $bamDir/*.sorted.bam; do
	fn_IDsorted=${f%.bam}.IDsorted.bam
	fn_trueBam="${fn_IDsorted%.bam}.true.bam"
	fn_map=$bedDir/`basename ${f%.sorted.bam}_mapping.csv`
	if [ ! -f $fn_map ]; then 
		echo $0 mapping file not found $fn_map
		exit 0
	fi
	fn_sorted=`echo $fn_trueBam | sed 's/.IDsorted//'`
	fn_out=$bamDir/logs/`basename ${f%.bam}`_sel_true.o

	if [ ! -f $fn_sorted ]; then 
		command="`dirname $0`/../eval_vcf/select_true_alignment $fn_IDsorted $fn_trueBam --map $fn_map; $samtools sort -l 9 $fn_trueBam ${fn_sorted%.bam}"
		echo "$command"
		submit "$command" --tag sel_true --work-dir `pwd` --log-dir $bamDir/logs/ --queue mpi01.q --fn-out $fn_out
	fi
done
wait_for_jobs sel_true


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

for fn_bam in `find $bamDir -name *.sorted.true.bam`; do
	if [ -s ${fn_bam%.bam}.karyotypic.bam ]; then 
		continue
	fi
	fn_out=$bamDir/logs/`basename ${fn_bam%.bam}`_sort_karyo.o
	command="source `find $SNV_DIR -name sort_bam.sh`; source `find $SNV_DIR -name paths.sh`; sort_bam $fn_bam $bamDir $statusDir --bam-ext .sorted.true.bam --human-karyotypic"
	submit "$command"  --tag sort_`basename ${fn_bam%.bam}` --log-dir $bamDir/logs/ --work-dir $SNV_DIR --queue mpi01.q --fn-out $fn_out
done
wait_for_jobs sort_


## important before merging:
## intersect each bam file with the corresponding bed file --> we exclude the losses and restrict the gains to the gain regions
# first link the entire bed file for the normal sample (the normal sample does not have losses or gains)
for g in 0 1; do
	if [[ ! -f $bedDir/${g}_NO_final.bed ]]; then
		ln -s $bed_file $bedDir/${g}_NO_final.bed
	fi
done

# now intersect each bam file with the corresponding region
for bam in $bamDir/*.sorted.true.karyotypic.bam; do
	matchedBed=$bedDir/`basename ${bam%.sorted.true.karyotypic.bam}`.bed
	if [[ ! -f $matchedBed ]]; then
		echo "Macthed bed file cannot be found: $bam"
		continue
	fi
	intersectedBam=${bam%.bam}.intersected.bam
	fn_out=${bamDir}/logs/intersect_`basename ${bam%.bam}`.o
	command="source `find $SNV_DIR -name paths.sh`; ${samtools_1_2} view -bh -L $matchedBed $bam > $intersectedBam "
	if [[ ! -f $intersectedBam ]]; then
		echo $command
		submit "$command"  --tag intersect_`basename ${bam%.bam}` --log-dir $bamDir/logs/ --work-dir $SNV_DIR --queue mpi01.q --fn-out $fn_out
	fi
done
wait_for_jobs intersect_
##


# echo merge the TU and NO bam files together
outBamTU=$bamDir/TU_final.true.bam
if [ -f $outBamTU ]; then
	echo "$outBamTU already exists."
else
	fn_out=$bamDir/logs/merge_job_`basename ${outBamTU%.bam}`.o
	command="$samtools merge $outBamTU $bamDir/*TU*.sorted.true.karyotypic.intersected.bam"
	submit "$command"  --tag merge_job_`basename ${outBamTU%.bam}` --log-dir $bamDir/logs/ --queue mpi01.q --fn-out $fn_out
fi

outBamNO=$bamDir/NO_final.true.bam
if [ -f $outBamNO ]; then
		echo "$outBamNO already exists."
else
	fn_out=$bamDir/logs/merge_job_`basename ${outBamNO%.bam}`.o
	command="$samtools merge $outBamNO $bamDir/*NO*.sorted.true.karyotypic.intersected.bam"
	submit "$command"  --tag merge_job_`basename ${outBamNO%.bam}` --log-dir $bamDir/logs/ --queue mpi01.q --fn-out $fn_out
fi
wait_for_jobs merge_job_


## echo replace read groups such that we only have one read group in the header (with only one sample name)
for curr_bam in $bamDir/NO_final.true.bam
do
	echo $curr_bam
	outputBam=${curr_bam%.true.bam}.RG.true.bam

	if [ -f $outputBam ]; then 
		continue
	fi

	RGID=`basename ${curr_bam%.bam}`
	RGLB=1
	RGPL="illumina"
	RGPU=1
	RGSM=`basename ${curr_bam%.bam} | tr '.' ' '  |  tr '_' ' ' | awk '{print $1}'`

	fn_out=$bamDir/logs/replace_RG`basename ${curr_bam%.bam}`.o
	command="source `find $SNV_DIR -name replaceReadGroups.sh`; source `find $SNV_DIR -name paths.sh`; replaceReadGroups $curr_bam $outputBam $RGID $RGLB $RGPL $RGPU $RGSM --check-fexist"
	submit "$command"  --tag replace_RG`basename ${curr_bam%.bam}` --log-dir $bamDir/logs/  --queue mpi01.q --fn-out $fn_out
done
wait_for_jobs replace_RG


echo create subsets of the bam files to merge with the TU file
bamFile=$bamDir/NO_final.true.bam
for fraction in 0.10 0.20 0.40 0.60; do 
	outBam=${bamFile%.true.bam}.`echo $fraction | cut -f 2 -d "."`perc.true.bam
	if [ -f $outBam ]; then 
		continue; 
	fi
	fn_out=$bamDir/logs/subset_trueNO_${fraction}.o
	command="source `find $SNV_DIR -name subsetSam.sh`; source `find $SNV_DIR -name paths.sh`; subsetSam $bamFile $outBam $fraction"
	submit "$command"  --tag subset1NO --log-dir $bamDir/logs/ --queue mpi01.q --fn-out $fn_out
done
wait_for_jobs subset1

# echo merge the TU final bam in the bamDir directory with the 10%, 20%, 40%, and 60% bam file
# to simulate normal contamination
for cont in 10 20 40 60
do
	outBamMerge=$bamDir/TU.wCont${cont}.final.true.bam
	if [ -f $outBamMerge ]; then
		echo "$outBamMerge already exists."
	else
		echo $bamFilesMerge
		fn_out=$bamDir/logs/merge_TU_`basename ${outBamMerge%.bam}`.o
		command="$samtools merge $outBamMerge $bamDir/TU_final.true.bam $bamDir/NO_final.${cont}perc.true.bam"
		submit "$command" --tag merge_TU_`basename ${outBamMerge%.bam}` --log-dir $bamDir/logs/ --queue mpi01.q --fn-out $fn_out
	fi
done
wait_for_jobs merge_TU_

## echo replace read groups such that we only have one read group in the header (with only one sample name)
for curr_bam in $bamDir/NO_final.true.bam $bamDir/TU.wCont[0-9]0.final.true.bam
do
	echo $curr_bam
	outputBam=${curr_bam%.true.bam}.RG.true.bam

	if [ -f $outputBam ]; then 
		continue
	fi

	RGID=`basename ${curr_bam%.bam}`
	RGLB=1
	RGPL="illumina"
	RGPU=1
	RGSM=`basename ${curr_bam%.bam} | tr '.' ' '  |  tr '_' ' ' | awk '{print $1}'`

	fn_out=$bamDir/logs/replace_RG`basename ${curr_bam%.bam}`.o
	command="source `find $SNV_DIR -name replaceReadGroups.sh`; source `find $SNV_DIR -name paths.sh`; replaceReadGroups $curr_bam $outputBam $RGID $RGLB $RGPL $RGPU $RGSM --check-fexist"
	submit "$command"  --tag replace_RG`basename ${curr_bam%.bam}` --log-dir $bamDir/logs/ --queue mpi01.q --fn-out $fn_out
done
wait_for_jobs replace_RG


for frac in 12 25 50 75 100; do

	all_cont=20
	if [ $frac -eq 50 ]; then 
		all_cont="10 20 40 60"
		echo run for all: $all_cont
	fi
	for cont in $all_cont; do
		bamFile=$bamDir/TU.wCont${cont}.final.RG.true.bam
		fraction=`echo "scale=6; $frac/1.${cont}/100" | bc`
		outBam=${bamFile%.true.bam}.${frac}perc.true.bam
		echo create file $outBam
		if [ -f $outBam ]; then 
			continue; 
		fi
		fn_out=$bamDir/logs/subset_for_`basename ${outBam%.bam}`.o
		command="source `find $SNV_DIR -name subsetSam.sh`; source `find $SNV_DIR -name paths.sh`; subsetSam $bamFile $outBam $fraction; $samtools index $outBam"
		submit "$command"  --tag subsetTU --log-dir $bamDir/logs/ --fn-out $fn_out
	done

	#Normal files
	bamFile=$bamDir/NO_final.RG.true.bam
	fraction=`echo "scale=6; $frac/100" | bc`
	outBam=${bamFile%.true.bam}.${frac}perc.true.bam
	if [ -f $outBam ]; then 
		continue; 
	fi
	fn_out=$bamDir/logs/subset_for_`basename ${outBam%.bam}`.o
	command="source `find $SNV_DIR -name subsetSam.sh`; source `find $SNV_DIR -name paths.sh`; subsetSam $bamFile $outBam $fraction; $samtools index $outBam"
	submit "$command"  --tag subsetNO --log-dir $bamDir/logs/  --fn-out $fn_out

done
#wait_for_jobs subset

