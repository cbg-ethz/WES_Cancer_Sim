#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`

for SAMVCF in `ls $base_dir/alignments_vsn_k20/variants_tuned/*sam*.gz | grep -v combined | grep -v NO | grep -v GERMLINE | grep -v eval | grep -v SOMATIC | grep -v Germline` `find $base_dir/alignments_vsn_k20/variants_default/ -name TU*samvar*.vcf.gz | grep -v eval | grep -v GERMLINE | grep -v SOMATIC | grep -v combined | grep -v SOMATIC_paired` `find $base_dir/alignments_sn_k1/variants_default/ -name TU*samvar*.vcf.gz | grep -v eval | grep -v GERMLINE | grep -v SOMATIC | grep -v combined | grep -v SOMATIC_paired`  `find $base_dir/alignments_vsn_k20/bam/variants_repeatedSubsampling/ -name TU*samvar*.vcf.gz | grep -v eval | grep -v GERMLINE | grep -v SOMATIC |  grep -v combined | grep -v SOMATIC_paired`; do
	NoVCF=${SAMVCF%_samvar_1_2.vcf.gz}_noE_samvar_1_2.vcf.NOsample.gz
	VCF_TUMORgermline=${SAMVCF%.vcf.gz}.GERMLINE.vcf
	VCF_TUMORsomatic=${SAMVCF%.vcf.gz}.SOMATIC.vcf
	if [ ! -f $VCF_TUMORgermline.gz ]; then
		command="java WriteDifference_TU_NO_asSomaticFilter $SAMVCF $NoVCF $VCF_TUMORgermline $VCF_TUMORsomatic; gzip $VCF_TUMORgermline ; gzip $VCF_TUMORsomatic"
		echo $command
		eval $command
	fi
done

# This is for both gatkUG and gatkHPCaller
for gatkVCF in `find $base_dir/alignments_vsn_k20/variants_default/ -name *gatk*.vcf | grep -v eval | grep -v GERMLINE | grep -v SOMATIC | grep -v combined | grep -v NO | grep -v paired ` `find $base_dir/alignments_sn_k1/variants_default/ -name *gatk*.vcf | grep -v eval | grep -v GERMLINE | grep -v SOMATIC | grep -v combined | grep -v NO | grep -v paired ` `find $base_dir/alignments_vsn_k20/bam/variants_repeatedSubsampling/ -name *gatk*.vcf | grep -v eval | grep -v GERMLINE | grep -v SOMATIC |  grep -v combined | grep -v NO | grep -v paired ` ; do
	NoVCF=`echo ${gatkVCF} | sed 's/.gatk/.NOsample.gatk/'`
	VCF_TUMORgermline=${gatkVCF%.vcf}.GERMLINE.vcf
	VCF_TUMORsomatic=${gatkVCF%.vcf}.SOMATIC.vcf
	if [ ! -f $VCF_TUMORgermline ]; then
		command="java WriteDifference_TU_NO_asSomaticFilter $gatkVCF $NoVCF $VCF_TUMORgermline $VCF_TUMORsomatic"
		echo $command
		eval $command
	else
		echo "$VCF_TUMORgermline already exists"
	fi
done

