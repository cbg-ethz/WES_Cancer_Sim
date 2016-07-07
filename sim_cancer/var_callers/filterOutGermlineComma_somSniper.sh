#!/bin/bash

## This needs to be done for somaticSniper
## somSniper:
currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`
for vcf in `find $base_dir/alignments_vsn_k20/variants_default/ -name *somaticSniper*.vcf | grep -v noComma | grep -v eval | grep -v GERMLINE | grep qual | grep -v combined | grep -v NO` `find $base_dir/alignments_sn_k1/variants_default/ -name *somaticSniper*.vcf | grep -v noComma | grep -v eval | grep -v GERMLINE | grep qual | grep -v combined | grep -v NO` `find $base_dir/alignments_vsn_k20/bam/variants_repeatedSubsampling/ -name *somaticSniper*.vcf | grep -v noComma | grep -v eval | grep -v GERMLINE | grep qual | grep -v combined | grep -v NO`  ; do
	outVCF=${vcf%.vcf}.noComma.vcf
	if [ ! -f $outVCF ]; then
		command="java Extract_Somatic_Mult_Alleles $vcf $outVCF"
		echo $command
		eval $command
	fi
done
exit 0 

