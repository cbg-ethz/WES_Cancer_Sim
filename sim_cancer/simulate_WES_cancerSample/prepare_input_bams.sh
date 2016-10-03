#!/bin/bash
currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`
source `find ${currScriptDir}/../ -name run_smmr.sh`

echo $gitDir
bamDir=$1

if [ -z "$bamDir" ]; then
	echo "Please provide bam directory"
	exit 0
fi

cd $bamDir

${samtools} view ccRCC_TU_mm.rmdup.bam -h | awk '{OFS="\t"; flag=$2; if (flag>=2048) flag-=2048; if (flag>=1024) flag-=1024; if (flag>=512) flag-=512; if (flag>=256) $2-=256;  print }' | ${samtools} view -bS '-' > ccRCC_TU_mm.ALL.rmdup.bam  
${samtools} index ccRCC_TU_mm.ALL.rmdup.bam

${samtools} view ccRCC_NO_mm.rmdup.bam -h | awk '{OFS="\t"; flag=$2; if (flag>=2048) flag-=2048; if (flag>=1024) flag-=1024; if (flag>=512) flag-=512; if (flag>=256) $2-=256;  print }' | ${samtools} view -bS '-' > ccRCC_NO_mm.ALL.rmdup.bam  
${samtools} index ccRCC_NO_mm.ALL.rmdup.bam 

run_smmr ccRCC_TU_mm.rmdup.bam
run_smmr ccRCC_NO_mm.rmdup.bam

mv ccRCC_NO_mm_smmr.rmdup.bai ccRCC_NO_mm_smmr.rmdup.bam.bai
mv ccRCC_TU_mm_smmr.rmdup.bai ccRCC_TU_mm_smmr.rmdup.bam.bai

cd -


