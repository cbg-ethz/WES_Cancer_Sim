#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`
source `find ${gitDir} -name paths.sh`
source `find $SNV_DIR -name utils.sh`
source `find $SNV_DIR -name submit.sh`

if [ -z $1 ]; then
	echo "Error: Please provide bam directory"
	exit 0
fi

bamDir=$1
outDir=$bamDir/repeatedSubsampling
mkdir -p $outDir
mkdir -p $outDir/logs

frac=50
cont=20
TUbamFile=$bamDir/TU.wCont${cont}.final.RG.bam
TUfraction=`echo "scale=6; $frac/1.${cont}/100" | bc`
echo "TUfraction=$TUfraction"
#Normal files
NObamFile=$bamDir/NO_final.RG.bam
NOfraction=`echo "scale=6; $frac/100" | bc`

myRandomSeeds=""
for repeat in `seq 1 10`; do
	myRandomSeeds="$myRandomSeeds $RANDOM"
done
echo $myRandomSeeds


for seed in $myRandomSeeds; do
	echo "Current seed = $seed"
	TUoutBam=$outDir/`basename ${TUbamFile%.bam}`.${frac}perc.${seed}.bam
	NOoutBam=$outDir/`basename ${NObamFile%.bam}`.${frac}perc.${seed}.bam
	echo create file $TUoutBam
	if [ -f $TUoutBam ]; then 
		continue; 
	fi
	command="source `find $SNV_DIR -name subsetSam.sh`; source `find $SNV_DIR -name paths.sh`; subsetSam $TUbamFile $TUoutBam $TUfraction --my-seed $seed; $samtools index $TUoutBam"
	submit "$command"  --tag subsetTU --queue mpi01.q --fn-out $outDir/logs/`basename ${TUoutBam%.bam}`.o

	if [ -f $NOoutBam ]; then 
		continue; 
	fi
	command="source `find $SNV_DIR -name subsetSam.sh`; source `find $SNV_DIR -name paths.sh`; subsetSam $NObamFile $NOoutBam $NOfraction --my-seed $seed; $samtools index $NOoutBam"
	submit "$command"  --tag subsetNO --queue mpi01.q --fn-out $outDir/logs/`basename ${NOoutBam%.bam}`.o
done


