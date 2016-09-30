#!/bin/bash

#percentFDR=10
#percentFDR=5
percentFDR=1

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`

SNVlist_percentFolder=$base_dir/alignments_vsn_k20/variants_default/`echo SNV_lists_fdr0${percentFDR}_${percentFDR}percent | sed 's/fdr010_/fdr10_/'`
outDir=${SNVlist_percentFolder}/overlaps/allOverlapsVennTop5ToolsExtractLists
varDir=$base_dir/alignments_vsn_k20/variants_default/

originalDeepSNVFile=TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf
originalSinvictFile=sinvict_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf
originalGatkUGFile=TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf
originalJointSNVMix2File=TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf
originalSomaticSniperFile=TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf


mkdir -p $outDir
cd $outDir
for file in gatkUG deep joint somSniper sinvict; do
	if [ ! -f ${file}.txt ]; then
		ln -s ${SNVlist_percentFolder}/overlaps/${file}.txt ${file}.txt
	fi
	if [ ! -f ${file}_pos.txt ]; then
		cat ${file}.txt | awk '{print $1,$2}' | tr ' ' '_' > ${file}_pos.txt
	fi
done
cd -


inputs=""
for file in gatkUG deep joint somSniper sinvict; do
	inputs="$inputs $outDir/${file}_pos.txt"
done

if [ ! -f $outDir/intersectAll.txt ]; then
	$currScriptDir/getSets.R $inputs
fi

for file in $outDir/*private*.txt; do
	firstTool=$( echo `basename $file` | tr '_' ' ' | awk '{print $1}')
	echo $firstTool
	secondTool=$( echo `basename $file` | tr '_' ' ' | awk '{print $3}')
	
	if [ "$secondTool" != "pos" -a "$secondTool" != "pos.txt" ]; then
		allTools="$firstTool $secondTool"
	else
		allTools="$firstTool"
	fi

	for tool in $allTools; do
		if [ "$tool"  == "deep" ]; then
			originalVCF=$varDir/${originalDeepSNVFile}
		elif [ "$tool" == "somSniper" ]; then
			originalVCF=$varDir/${originalSomaticSniperFile}
		elif [ "$tool" == "gatkUG" ]; then
			originalVCF=$varDir/${originalGatkUGFile}
		elif [ "$tool" == "sinvict" ]; then
			originalVCF=$varDir/${originalSinvictFile}
		elif [ "$tool" == "joint" ]; then
			originalVCF=$varDir/${originalJointSNVMix2File}
		else
			echo "Error: tool=$tool"
			exit 0
		fi
		outVCF=${file%.txt}_original_${tool}.vcf
		if [ ! -f $outVCF ]; then
			command="java createOriginalVCF $file $originalVCF $outVCF"
			echo $command
			eval $command
		fi
	done
done

for file in $outDir/intersectAll.txt $outDir/joint_overlap_sinvict_private_pos.txt; do
	originalVCF=$varDir/${originalJointSNVMix2File}
	outVCF=${file%.txt}_original_joint.vcf
	if [ ! -f $outVCF ]; then
		command="java createOriginalVCF $file $originalVCF $outVCF"
		echo $command
		eval $command
	fi
done

for file in ${SNVlist_percentFolder}/*.vcf; do
	cat $file | grep -v ^# | awk '{print $1"_"$2}' | sort | uniq > ${file}_pos
	tool=`basename $file`
	if [[ "$tool"  == *deep* ]]; then
		originalVCF=$varDir/${originalDeepSNVFile}
	elif [[ "$tool" == *somaticSniper* ]]; then
		originalVCF=$varDir/${originalSomaticSniperFile}
	elif [[ "$tool" == *gatk_SNVs* ]]; then
		originalVCF=$varDir/${originalGatkUGFile}
	elif [[ "$tool" == *sinvict* ]]; then
		originalVCF=$varDir/${originalSinvictFile}
	elif [[ "$tool" == *joint* ]]; then
		originalVCF=$varDir/${originalJointSNVMix2File}
	else
		echo "Error: tool=$tool"
		#exit 0
		continue
	fi
	outVCF=${file%.vcf}_original_VCF.vcf
	if [ ! -f $outVCF ]; then
		command="java createOriginalVCF ${file}_pos $originalVCF $outVCF"
		echo $command
		eval $command
	fi
done

