#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`

for somaticSniperVCF in `find $base_dir/alignments_*/variants_default/ -name *somaticSniper_SNVs_Raw.vcf` `find $base_dir/alignments_vsn_k20/bam/variants_repeatedSubsampling -name *somaticSniper_SNVs_Raw.vcf`
do
	if [ ! -f ${somaticSniperVCF%.vcf}_qual.vcf ]; then
		echo "$somaticSniperVCF"

		command="java RewriteSomaticSniperQualField $somaticSniperVCF > ${somaticSniperVCF%.vcf}_qual.vcf "
		echo $command
		eval $command
	fi
done



for varscanVCF in `find $base_dir/alignments_*/variants*/ -name *varscan2.txt.snp.Somatic.vcf` `find $base_dir/alignments_vsn_k20/bam/variants_repeatedSubsampling -name *varscan2.txt.snp.Somatic.vcf`
do
	if [ ! -f ${varscanVCF%.vcf}_qual.vcf ]; then
		echo "VarScan2"
		command="java -jar $currScriptDir/RewriteVarScanVCF.jar $varscanVCF > ${varscanVCF%.vcf}_qual.vcf"
		echo $command
		eval $command
	fi
done
