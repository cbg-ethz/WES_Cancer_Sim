#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`

## rewrite/filter the new variant calls:
# somaticsniper and varscan: quality field
for somaticSniperVCF in `find $base_dir_wCNVs/alignments_*/variants_default/ -name *somaticSniper_SNVs_Raw.vcf`
do
	if [ ! -f ${somaticSniperVCF%.vcf}_qual.vcf ]; then
		echo "$somaticSniperVCF"

		command="java RewriteSomaticSniperQualField $somaticSniperVCF > ${somaticSniperVCF%.vcf}_qual.vcf "
		echo $command
		eval $command
	fi
done
for varscanVCF in `find $base_dir_wCNVs/alignments_*/variants*/ -name *varscan2.txt.snp.Somatic.vcf`
do
	if [ ! -f ${varscanVCF%.vcf}_qual.vcf ]; then
		echo "VarScan2"
		command="java -jar RewriteVarScanVCF.jar $varscanVCF > ${varscanVCF%.vcf}_qual.vcf"
		echo $command
		eval $command
	fi
done


## gatk and samtools: filter out germline:
for SAMVCF in `find $base_dir_wCNVs/alignments_vsn_k20/variants_default/ -name TU*samvar*.vcf.gz | grep -v eval_1000 | grep -v GERMLINE | grep -v SOMATIC | grep -v combined | grep -v SOMATIC_paired`; do
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
for gatkVCF in `find $base_dir_wCNVs/alignments_vsn_k20/variants_default/ -name *gatk*.vcf | grep -v eval_1000 | grep -v GERMLINE | grep -v SOMATIC | grep -v combined | grep -v NO | grep -v paired `; do
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

# format HP Caller:
for vcf in `find $base_dir_wCNVs/alignments_*/variants_*/ -name *HPCaller*.vcf | grep -v rewritten | grep SOMATIC | grep -v paired`
do
        if [ ! -f ${vcf%.vcf}.rewritten.vcf ]; then
                command="java FormatHPCallerGATK_TUonly ${vcf} > ${vcf%.vcf}.rewritten.vcf"
                echo $command
                eval $command
        fi
done
# som sniper comma
for vcf in `find $base_dir_wCNVs/alignments_vsn_k20/variants_default/ -name *somaticSniper*.vcf | grep -v noComma | grep -v eval_1000 | grep -v GERMLINE | grep qual | grep -v combined | grep -v NO`; do
        outVCF=${vcf%.vcf}.noComma.vcf
        if [ ! -f $outVCF ]; then
                command="java Extract_Somatic_Mult_Alleles $vcf $outVCF"
                echo $command
                eval $command
        fi
done


