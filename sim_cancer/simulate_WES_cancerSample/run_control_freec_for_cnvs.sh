#!/bin/bash
currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`
source `find ${currScriptDir}/../ -name run_smmr.sh`
source `find $SNV_DIR -name paths.sh`
source `find $SNV_DIR -name utils.sh`
source `find $SNV_DIR -name submit.sh`
source `find $SNV_DIR -name excavator_cnvs.sh`

echo $gitDir
bamDir=$1
outDir=$2

if [ -z "$bamDir" -o -z "$outDir" ]; then
	echo "Please provide bam directory and out dir"
	exit 0
fi

tumorBam=$bamDir/19_4845-M_TU1_C3_mm_smmr.rmdup.bam
normalBam=$bamDir/21_4845-M_NO_E3_mm_smmr.rmdup.bam

mkdir -p $outDir


fn_config=$outDir/Control-FREEC.config

echo "[general]:" > $fn_config
echo "BedGraphOutput=TRUE" >> $fn_config
echo "bedtools=$bedtools"
echo "chrLenFile = $chrLen" >> $fn_config
echo "contaminationAdjustment=TRUE" >> $fn_config
echo "maxThreads=4" >> $fn_config
echo "outputDir = $outDir" >> $fn_config
echo "ploidy=2" >> $fn_config
echo "sex=XX" >> $fn_config
echo "window=0" >> $fn_config
echo "printNA=FALSE" >> $fn_config
echo "breakPointType=4" >> $fn_config
echo "chrFiles=$pathToChrFiles" >> $fn_config
echo "breakPointThreshold=1.2" >> $fn_config
echo "noisyData=TRUE" >> $fn_config
echo "readCountThreshold=50" >> $fn_config
echo "gemMappabilityFile=$gemMapFile" >> $fn_config
echo "minCNAlength=50" >> $fn_config
echo "uniqueMatch=TRUE" >> $fn_config

echo "[sample]" >> $fn_config
echo " " >> $fn_config
echo "mateFile = $tumorBam" >> $fn_config
echo "inputFormat=bam" >> $fn_config
echo "mateOrientation=FR" >> $fn_config
echo " " >> $fn_config
echo "[control]" >> $fn_config
echo " " >> $fn_config
echo "mateFile = $normalBam" >> $fn_config
echo "inputFormat=bam" >> $fn_config
echo "mateOrientation=FR" >> $fn_config
echo " " >> $fn_config
echo "[target]" >> $fn_config
echo " " >> $fn_config
echo "captureRegions = $bed_file" >> $fn_config

command="$freec -conf $fn_config"
echo $command
eval $command
