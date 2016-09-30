#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`
source `find ${gitDir} -name paths.sh`
tools_dir=${gitDir}
source `find ${tools_dir} -name genome.sh`



cnvListDir=$1
variantListDir=$2

if [ -z "$1" -o -z "$2" ]; then
	echo "cnvListDir and variantListDir needs to be provided"
	exit 0
fi


# define outDir
outDir=$cnvListDir/intersectedVariants/
mkdir -p $outDir

# for the losses, take all variants which are not in the loss regions:
# first for the copy number 0:
	## no loss:
	nodeInTree=7891011121314_0_TU
	ln -s $variantListDir/overlap_7891011121314_0_TU.txt $outDir/node_${nodeInTree}.txt

	## loss file: cnv_losses_7_0_TU_8_0_TU_9_0_TU_10_0_TU.bed
	lossFile=$cnvListDir/cnv_losses_7_0_TU_8_0_TU_9_0_TU_10_0_TU.bed
	formerNodeInTree=7891011121314_0_TU
	nodeInTree=7_0_TU_8_0_TU_9_0_TU_10_0_TU
	# file with variants for that node: private_7891011121314_0_TU.txt
	$intersectBed -a $variantListDir/private_${formerNodeInTree}.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	lossFile=$cnvListDir/cnv_losses_7_0_TU_8_0_TU_9_0_TU_10_0_TU.bed
	nodeInTree=7_0_TU_8_0_TU
	$intersectBed -a $variantListDir/private_7_0_TU_8_0_TU_9_0_TU_10_0_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	nodeInTree=9_0_TU_10_0_TU
	$intersectBed -a $variantListDir/private_9_0_TU_10_0_TU_7_0_TU_8_0_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	nodeInTree=8_0_TU
	$intersectBed -a $variantListDir/private_8_0_TU_minus_7_0_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	nodeInTree=9_0_TU
	$intersectBed -a $variantListDir/private_9_0_TU_minus_10_0_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	
	nodeInTree=10_0_TU
	$intersectBed -a $variantListDir/private_10_0_TU_minus_9_0_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	## additional loss file: cnv_losses_7_0_TU.bed
	formerNodeInTree=7_0_TU_8_0_TU
	nodeInTree=7_0_TU
	# file with variants for that node: private_7_0_TU_minus_8_0_TU.txt
	# inherit the loss regions which are deleted:
	combinedLossFile=$cnvListDir/cnv_losses_inclInherited_${nodeInTree}.bed
	cat $cnvListDir/cnv_losses_7_0_TU_8_0_TU_9_0_TU_10_0_TU.bed $cnvListDir/cnv_losses_${nodeInTree}.bed > $combinedLossFile
	$intersectBed -a $variantListDir/private_7_0_TU_minus_8_0_TU.txt -b $combinedLossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	## no loss:
	nodeInTree=11_0_TU_12_0_TU_13_0_TU_14_0_TU
	ln -s $variantListDir/private_1112131478910_0_TU.txt $outDir/node_${nodeInTree}.txt	

	## loss file: cnv_losses_11_0_TU_12_0_TU.bed
	lossFile=$cnvListDir/cnv_losses_11_0_TU_12_0_TU.bed
	formerNodeInTree=11_0_TU_12_0_TU_13_0_TU_14_0_TU
	nodeInTree=11_0_TU_12_0_TU
	# file with variants for that node: private_11_0_TU_12_0_TU_13_0_TU_14_0_TU.txt
	$intersectBed -a $variantListDir/private_${formerNodeInTree}.txt -b $cnvListDir/cnv_losses_${nodeInTree}.bed -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	## no loss:
	nodeInTree=13_0_TU_14_0_TU
	ln -s $variantListDir/private_13_0_TU_14_0_TU_11_0_TU_12_0_TU.txt $outDir/node_${nodeInTree}.txt

	lossFile=$cnvListDir/cnv_losses_11_0_TU_12_0_TU.bed

	nodeInTree=11_0_TU
	$intersectBed -a $variantListDir/private_11_0_TU_minus_12_0_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	nodeInTree=12_0_TU
	$intersectBed -a $variantListDir/private_12_0_TU_minus_11_0_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt	

	## no loss:
	nodeInTree=13_0_TU
	ln -s $variantListDir/private_13_0_TU_minus_14_0_TU.txt $outDir/node_${nodeInTree}.txt
	
	## loss file: cnv_losses_14_0_TU.bed
	lossFile=$cnvListDir/cnv_losses_14_0_TU.bed
	formerNodeInTree=13_0_TU_14_0_TU
	nodeInTree=14_0_TU
	# file with variants for that node: private_14_0_TU_minus_13_0_TU.txt
	$intersectBed -a $variantListDir/private_14_0_TU_minus_13_0_TU.txt -b $cnvListDir/cnv_losses_${nodeInTree}.bed -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

# then for the copy number 1:
	nodeInTree=7891011121314_1_TU
	lossFile=$cnvListDir/cnv_losses_7891011121314_1_TU.bed
	$intersectBed -a $variantListDir/overlap_${nodeInTree}.txt -b $cnvListDir/cnv_losses_${nodeInTree}.bed -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	formerNodeInTree=7891011121314_1_TU
	nodeInTree=7_1_TU_8_1_TU_9_1_TU_10_1_TU 
	# file with variants for that node: private_7891011121314_1_TU.txt
	# additional loss file: $cnvListDir/cnv_losses_7_1_TU_8_1_TU_9_1_TU_10_1_TU.bed
	cat $cnvListDir/cnv_losses_${formerNodeInTree}.bed $cnvListDir/cnv_losses_${nodeInTree}.bed > $cnvListDir/cnv_losses_inclInherited_${nodeInTree}.bed
	$intersectBed -a $variantListDir/private_${formerNodeInTree}.txt -b $cnvListDir/cnv_losses_inclInherited_${nodeInTree}.bed -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	lossFile=$cnvListDir/cnv_losses_inclInherited_7_1_TU_8_1_TU_9_1_TU_10_1_TU.bed

	nodeInTree=7_1_TU_8_1_TU
	$intersectBed -a $variantListDir/private_7_1_TU_8_1_TU_9_1_TU_10_1_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	nodeInTree=9_1_TU_10_1_TU
	$intersectBed -a $variantListDir/private_9_1_TU_10_1_TU_7_1_TU_8_1_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	formerNodeInTree=7_1_TU_8_1_TU
	nodeInTree=7_1_TU
	combinedLossFile=$cnvListDir/cnv_losses_inclInherited_${nodeInTree}.bed
	cat $cnvListDir/cnv_losses_inclInherited_7_1_TU_8_1_TU_9_1_TU_10_1_TU.bed $cnvListDir/cnv_losses_${nodeInTree}.bed > $combinedLossFile
	# file with variants for that node: private_7_1_TU_minus_8_1_TU.txt
	$intersectBed -a $variantListDir/private_7_1_TU_minus_8_1_TU.txt -b $combinedLossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	lossFile=$cnvListDir/cnv_losses_inclInherited_7_1_TU_8_1_TU_9_1_TU_10_1_TU.bed

	nodeInTree=8_1_TU
	$intersectBed -a $variantListDir/private_8_1_TU_minus_7_1_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt
	
	nodeInTree=9_1_TU
	$intersectBed -a $variantListDir/private_9_1_TU_minus_10_1_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt	

	nodeInTree=10_1_TU
	$intersectBed -a $variantListDir/private_10_1_TU_minus_9_1_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	formerNodeInTree=7891011121314_1_TU
	nodeInTree=11_1_TU_12_1_TU_13_1_TU_14_1_TU
	# file with variants for that node: private_1112131478910_1_TU.txt
	combinedLossFile=$cnvListDir/cnv_losses_inclInherited_${nodeInTree}.bed
	cat $cnvListDir/cnv_losses_${formerNodeInTree}.bed $cnvListDir/cnv_losses_${nodeInTree}.bed > $combinedLossFile
	$intersectBed -a $variantListDir/private_1112131478910_1_TU.txt -b $combinedLossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	lossFile=$cnvListDir/cnv_losses_inclInherited_11_1_TU_12_1_TU_13_1_TU_14_1_TU.bed
	
	nodeInTree=11_1_TU_12_1_TU
	$intersectBed -a $variantListDir/private_11_1_TU_12_1_TU_13_1_TU_14_1_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	nodeInTree=13_1_TU_14_1_TU
	$intersectBed -a $variantListDir/private_13_1_TU_14_1_TU_11_1_TU_12_1_TU.txt  -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	nodeInTree=11_1_TU
	$intersectBed -a $variantListDir/private_11_1_TU_minus_12_1_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt
	
	nodeInTree=12_1_TU
	$intersectBed -a $variantListDir/private_12_1_TU_minus_11_1_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	nodeInTree=13_1_TU
	$intersectBed -a $variantListDir/private_13_1_TU_minus_14_1_TU.txt -b $lossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt

	formerNodeInTree=13_1_TU_14_1_TU
	nodeInTree=14_1_TU
	# file with variants for that node: private_14_1_TU_minus_13_1_TU.txt
	combinedLossFile=$cnvListDir/cnv_losses_inclInherited_${nodeInTree}.bed
	cat $cnvListDir/cnv_losses_inclInherited_11_1_TU_12_1_TU_13_1_TU_14_1_TU.bed $cnvListDir/cnv_losses_${nodeInTree}.bed > $combinedLossFile
	$intersectBed -a $variantListDir/private_14_1_TU_minus_13_1_TU.txt -b $combinedLossFile -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt



