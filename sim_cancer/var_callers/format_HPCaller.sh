#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`

for vcf in `ls $base_dir/alignments_vsn_k20/bam/variants_repeatedSubsampling/*gatkHPCaller*.vcf | grep SOMATIC | grep -v paired | grep -v rewritten` `find $base_dir/alignments_*/variants_*/ -name *HPCaller*.vcf | grep -v rewritten | grep SOMATIC | grep -v paired`
do
	if [ ! -f ${vcf%.vcf}.rewritten.vcf ]; then
		command="java FormatHPCallerGATK_TUonly ${vcf} > ${vcf%.vcf}.rewritten.vcf"
		echo $command
		eval $command
	fi
done
