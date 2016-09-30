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
outDir=$variantListDir/inclGains/
mkdir -p $outDir
## for each subclone that inherits CNVs: decide randomly, which mutations stay on the "main chromosomes", and which mutations are placed only on the CNV gain regions



# extract the region of gain from the node where it happens - all mutations are "duplicated" - including the ones that are in ancestor clones!!
nodeInTree=7_0_TU_8_0_TU
gainFile=$cnvListDir/cnv_gains_7_0_TU_8_0_TU.bed
$intersectBed -a $variantListDir/node_${nodeInTree}_withoutLossRegions.txt -b $gainFile > $outDir/node_${nodeInTree}_gain.txt ## Note: no mutations in this region!

# including the ones that are in ancestor clones!
AncestorNodeInTree=7_0_TU_8_0_TU_9_0_TU_10_0_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
AncestorNodeInTree=7891011121314_0_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt

# inherit this gain and distribute new mutations within the gain region randomly onto the "normal" genome or the gain region
nodeInTree=7_0_TU
gainFile=$cnvListDir/cnv_gains_7_0_TU_8_0_TU.bed
percentage=0.5
## for each subclone that inherits CNVs: decide randomly, which mutations stay on the "main chromosomes", and which mutations are placed only on the CNV gain regions
$intersectBed -a $variantListDir/node_${nodeInTree}_withoutLossRegions.txt -b $gainFile | awk "BEGIN {srand()} { if (rand() <= $percentage) print }"  > $outDir/node_${nodeInTree}_gain.txt ## Note: no mutations in this region!

# including the ones that are in ancestor clones!
AncestorNodeInTree=7_0_TU_8_0_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
AncestorNodeInTree=7_0_TU_8_0_TU_9_0_TU_10_0_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
AncestorNodeInTree=7891011121314_0_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt


# the rest of the mutations, which did not get put into the gain will be written into the mutation list for the "normal" genome:
if [ `cat $outDir/node_${nodeInTree}_gain.txt | wc -l ` == 0 ]; then # if it is empty, we link the whole variant file (no mutations in gain region)
	ln -s $variantListDir/node_${nodeInTree}_withoutLossRegions.txt $outDir/node_${nodeInTree}.txt
else # if there are variants in there
	java SNV_Overlapper_Ordered_PlusAlleles $variantListDir/node_${nodeInTree}_withoutLossRegions.txt $outDir/node_${nodeInTree}_gain.txt $outDir/node_${nodeInTree}_inGainAndOrigMainFile.txt  $overlapFile1File2 $outDir/node_${nodeInTree}.txt
fi


# inherit this gain and distribute new mutations within the gain region randomly onto the "normal" genome or the gain region
nodeInTree=8_0_TU
gainFile=$cnvListDir/cnv_gains_7_0_TU_8_0_TU.bed
percentage=0.5
$intersectBed -a $variantListDir/node_${nodeInTree}_withoutLossRegions.txt -b $gainFile | awk "BEGIN {srand()} { if (rand() <= $percentage) print }"  > $outDir/node_${nodeInTree}_gain.txt #Note: no mutations in this region!

# including the ones that are in ancestor clones!
AncestorNodeInTree=7_0_TU_8_0_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
AncestorNodeInTree=7_0_TU_8_0_TU_9_0_TU_10_0_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
AncestorNodeInTree=7891011121314_0_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt

# the rest of the mutations, which did not get put into the gain will be written into the mutation list for the "normal" genome:
if [ `cat $outDir/node_${nodeInTree}_gain.txt | wc -l ` == 0 ]; then # if it is empty, we link the whole variant file (no mutations in gain region)
	ln -s $variantListDir/node_${nodeInTree}_withoutLossRegions.txt $outDir/node_${nodeInTree}.txt
else # if there are variants in there
	java SNV_Overlapper_Ordered_PlusAlleles $variantListDir/node_${nodeInTree}_withoutLossRegions.txt $outDir/node_${nodeInTree}_gain.txt $outDir/node_${nodeInTree}_inGainAndOrigMainFile.txt  $overlapFile1File2 $outDir/node_${nodeInTree}.txt
fi

## Now the other branch of the tree:
## here, it happened as a homozygous event:
for copy in 0 1; do
	# extract the region of gain from the node where it happens - all mutations are "duplicated"
	nodeInTree=11_${copy}_TU_12_${copy}_TU
	gainFile=$cnvListDir/cnv_gains_11_${copy}_TU_12_${copy}_TU.bed
	$intersectBed -a $variantListDir/node_${nodeInTree}_withoutLossRegions.txt -b $gainFile > $outDir/node_${nodeInTree}_gain.txt

	# including the ones that are in ancestor clones!
	AncestorNodeInTree=11_${copy}_TU_12_${copy}_TU_13_${copy}_TU_14_${copy}_TU
	if [ -f $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt ]; then
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	else
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	fi
	AncestorNodeInTree=7891011121314_${copy}_TU
	if [ -f $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt ]; then
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	else
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	fi

	# inherit this gain and distribute new mutations within the gain region randomly onto the "normal" genome or the gain region
	nodeInTree=11_${copy}_TU
	gainFile=$cnvListDir/cnv_gains_11_${copy}_TU_12_${copy}_TU.bed
	percentage=0.5
	$intersectBed -a $variantListDir/node_${nodeInTree}_withoutLossRegions.txt -b $gainFile | awk "BEGIN {srand()} { if (rand() <= $percentage) print }"  > $outDir/node_${nodeInTree}_gain.txt 
	# including the ones that are in ancestor clones!
	AncestorNodeInTree=11_${copy}_TU_12_${copy}_TU
	if [ -f $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt ]; then
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	else
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	fi
	AncestorNodeInTree=11_${copy}_TU_12_${copy}_TU_13_${copy}_TU_14_${copy}_TU
	if [ -f $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt ]; then
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	else
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	fi
	AncestorNodeInTree=7891011121314_${copy}_TU
	if [ -f $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt ]; then
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	else
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	fi

	# the rest of the mutations, which did not get put into the gain will be written into the mutation list for the "normal" genome:
	if [ `cat $outDir/node_${nodeInTree}_gain.txt | wc -l ` == 0 ]; then # if it is empty, we link the whole variant file (no mutations in gain region)
		ln -s $variantListDir/node_${nodeInTree}_withoutLossRegions.txt $outDir/node_${nodeInTree}.txt
	else # if there are variants in there
		java SNV_Overlapper_Ordered_PlusAlleles $variantListDir/node_${nodeInTree}_withoutLossRegions.txt $outDir/node_${nodeInTree}_gain.txt $outDir/node_${nodeInTree}_inGainAndOrigMainFile.txt  $overlapFile1File2 $outDir/node_${nodeInTree}.txt
	fi
	# additionally, 11_${copy}_TU has another gain, which happens after the mutations:
	additionalGainFile=$cnvListDir/cnv_gains_11_${copy}_TU.bed
	# append all mutations which are in this region to the gain file
	$intersectBed -a $variantListDir/node_${nodeInTree}_withoutLossRegions.txt -b $additionalGainFile >> $outDir/node_${nodeInTree}_gain.txt


	# inherit this gain and distribute new mutations within the gain region randomly onto the "normal" genome or the gain region
	nodeInTree=12_${copy}_TU
	gainFile=$cnvListDir/cnv_gains_11_${copy}_TU_12_${copy}_TU.bed
	percentage=0.5
	$intersectBed -a $variantListDir/node_${nodeInTree}_withoutLossRegions.txt -b $gainFile | awk "BEGIN {srand()} { if (rand() <= $percentage) print }"  > $outDir/node_${nodeInTree}_gain.txt
	# including the ones that are in ancestor clones!
	AncestorNodeInTree=11_${copy}_TU_12_${copy}_TU
	if [ -f $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt ]; then
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	else
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	fi
	AncestorNodeInTree=11_${copy}_TU_12_${copy}_TU_13_${copy}_TU_14_${copy}_TU
	if [ -f $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt ]; then
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	else
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	fi
	AncestorNodeInTree=7891011121314_${copy}_TU
	if [ -f $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt ]; then
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	else
		$intersectBed -a $variantListDir/node_${AncestorNodeInTree}.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
	fi

	# the rest of the mutations, which did not get put into the gain will be written into the mutation list for the "normal" genome:
	if [ `cat $outDir/node_${nodeInTree}_gain.txt | wc -l ` == 0 ]; then # if it is empty, we link the whole variant file (no mutations in gain region)
		ln -s $variantListDir/node_${nodeInTree}_withoutLossRegions.txt $outDir/node_${nodeInTree}.txt
	else
		java SNV_Overlapper_Ordered_PlusAlleles $variantListDir/node_${nodeInTree}_withoutLossRegions.txt $outDir/node_${nodeInTree}_gain.txt $outDir/node_${nodeInTree}_inGainAndOrigMainFile.txt  $overlapFile1File2 $outDir/node_${nodeInTree}.txt
	fi
done


## whole chromosome gain:
# all of chromosome 17:
echo -e "chr17\t1\t81195210"

# extract the region of gain from the node where it happens - all mutations are "duplicated"
nodeInTree=13_1_TU_14_1_TU
gainFile=$cnvListDir/cnv_gains_13_1_TU_14_1_TU_allCHR17.bed
echo -e "chr17\t1\t81195210" > $gainFile
$intersectBed -a $variantListDir/node_${nodeInTree}_withoutLossRegions.txt -b $gainFile > $outDir/node_${nodeInTree}_gain.txt

# including the ones that are in ancestor clones!
AncestorNodeInTree=11_1_TU_12_1_TU_13_1_TU_14_1_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
AncestorNodeInTree=7891011121314_1_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt


# inherit this gain and distribute new mutations within the gain region randomly onto the "normal" genome or the gain region
nodeInTree=13_1_TU
gainFile=$cnvListDir/cnv_gains_13_1_TU_14_1_TU_allCHR17.bed
percentage=0.5
## for each subclone that inherits CNVs: decide randomly, which mutations stay on the "main chromosomes", and which mutations are placed only on the CNV gain regions
$intersectBed -a $variantListDir/node_${nodeInTree}_withoutLossRegions.txt -b $gainFile | awk "BEGIN {srand()} { if (rand() <= $percentage) print }"  > $outDir/node_${nodeInTree}_gain.txt ## Note: no mutations in this region!

# including the ones that are in ancestor clones!
AncestorNodeInTree=13_1_TU_14_1_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
AncestorNodeInTree=11_1_TU_12_1_TU_13_1_TU_14_1_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
AncestorNodeInTree=7891011121314_1_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt

# the rest of the mutations, which did not get put into the gain will be written into the mutation list for the "normal" genome:
if [ `cat $outDir/node_${nodeInTree}_gain.txt | wc -l ` == 0 ]; then # if it is empty, we link the whole variant file (no mutations in gain region)
	ln -s $variantListDir/node_${nodeInTree}_withoutLossRegions.txt $outDir/node_${nodeInTree}.txt
else # if there are variants in there
	java SNV_Overlapper_Ordered_PlusAlleles $variantListDir/node_${nodeInTree}_withoutLossRegions.txt $outDir/node_${nodeInTree}_gain.txt $outDir/node_${nodeInTree}_inGainAndOrigMainFile.txt  $overlapFile1File2 $outDir/node_${nodeInTree}.txt
fi

# inherit this gain and distribute new mutations within the gain region randomly onto the "normal" genome or the gain region
nodeInTree=14_1_TU
gainFile=$cnvListDir/cnv_gains_13_1_TU_14_1_TU_allCHR17.bed
percentage=0.5
$intersectBed -a $variantListDir/node_${nodeInTree}_withoutLossRegions.txt -b $gainFile | awk "BEGIN {srand()} { if (rand() <= $percentage) print }"  > $outDir/node_${nodeInTree}_gain.txt #Note: no mutations in this region!

# including the ones that are in ancestor clones!
AncestorNodeInTree=13_1_TU_14_1_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
AncestorNodeInTree=11_1_TU_12_1_TU_13_1_TU_14_1_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt
AncestorNodeInTree=7891011121314_1_TU
$intersectBed -a $variantListDir/node_${AncestorNodeInTree}_withoutLossRegions.txt -b $gainFile >> $outDir/node_${nodeInTree}_gain.txt

# the rest of the mutations, which did not get put into the gain will be written into the mutation list for the "normal" genome:
if [ `cat $outDir/node_${nodeInTree}_gain.txt | wc -l ` == 0 ]; then # if it is empty, we link the whole variant file (no mutations in gain region)
	ln -s $variantListDir/node_${nodeInTree}_withoutLossRegions.txt $outDir/node_${nodeInTree}.txt
else # if there are variants in there
	java SNV_Overlapper_Ordered_PlusAlleles $variantListDir/node_${nodeInTree}_withoutLossRegions.txt $outDir/node_${nodeInTree}_gain.txt $outDir/node_${nodeInTree}_inGainAndOrigMainFile.txt  $overlapFile1File2 $outDir/node_${nodeInTree}.txt
fi


