#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`
source `find ${gitDir} -name paths.sh`
tools_dir=${gitDir}
source `find ${tools_dir} -name genome.sh`


cnvListDir=$1
variantListDir=$2

if [ -z "$1" -o -z "$2" ]; then
	echo "cnvListDir, variantListDir needs to be provided"
	exit 0
fi

# define outDir
outDir=$variantListDir/finalListsAndRegions/
mkdir -p $outDir

# define vcf header:
vcfHeader="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample"

for copy in 0 1; do
	for ThiSubclone in 7 8 9 10 11 12 13 14; do
		subclone=${ThiSubclone}_${copy}_TU
		echo $subclone
		finalMutFileUnordered=$outDir/${ThiSubclone}_${copy}_TU_TEMP.txt
		finalMutFileIntersected=$outDir/${ThiSubclone}_${copy}_TU_TEMP_intersected.txt
		finalMutFile=$outDir/${ThiSubclone}_${copy}_TU.vcf

		# combine and sort
		if [ "$ThiSubclone" == "14" ]; then
			cat `ls $variantListDir/*${ThiSubclone}_${copy}_TU*.txt | grep -v gain` | sort -k1,1 -k2n > $finalMutFileUnordered
		else
			cat $variantListDir/node_7891011121314_${copy}_TU.txt `ls $variantListDir/*${ThiSubclone}_${copy}_TU*.txt | grep -v gain` | sort -k1,1 -k2n > $finalMutFileUnordered
		fi

		echo "make bed file"
		# determine bed file which contains all losses that this subclone has
		allLossRegionsBed=$outDir/${ThiSubclone}_${copy}_TU.bed
		if [ "$subclone" != "13_0_TU" ]; then
			cat $cnvListDir/cnv_losses*${ThiSubclone}_${copy}_TU* | sort -k1,1 -k2n | uniq > $outDir/${ThiSubclone}_${copy}_TU_TEMP.bed
			if [ -f $cnvListDir/cnv_losses_7891011121314_${copy}_TU.bed ]; then
				cat $cnvListDir/cnv_losses_7891011121314_${copy}_TU.bed >> $outDir/${ThiSubclone}_${copy}_TU_TEMP.bed
			fi
			cat $outDir/${ThiSubclone}_${copy}_TU_TEMP.bed | sort -k1,1 -k2n | uniq > $outDir/${ThiSubclone}_${copy}_TU_TEMP.bed.sortedAndUnique.bed
		else
			touch $outDir/${ThiSubclone}_${copy}_TU_TEMP.bed.sortedAndUnique.bed
		fi
		
		# create final bed file as intersection between the targeted regions and the regions which are NOT (-v) lost 
		if [ -s $outDir/${ThiSubclone}_${copy}_TU_TEMP.bed.sortedAndUnique.bed ]; then
			$intersectBed -a $bed_file -b $outDir/${ThiSubclone}_${copy}_TU_TEMP.bed.sortedAndUnique.bed -v > $allLossRegionsBed

			## regions which are only lost further down in the tree: all mutations from ancestor nodes need to be removed as well in these regions when adding up all mutations for a leaf --> therefore:
			# finally remove all mutations which were inherited but are in regions that were lost later:
			$intersectBed -a $finalMutFileUnordered -b $outDir/${ThiSubclone}_${copy}_TU_TEMP.bed.sortedAndUnique.bed -v > $finalMutFileIntersected

		else # empty file in the case of 13_0_TU
			ln -s $bed_file $allLossRegionsBed
			ln -s $finalMutFileUnordered $finalMutFileIntersected
		fi

		# add header and mutations to create final vcf file
		echo -e $vcfHeader > $finalMutFile
		cat $finalMutFileIntersected >> $finalMutFile
	done
done

# now for the gain files: the leaf subclone gain files already contain all mutations from the ancestor clones!
for subclone in 7_0_TU 8_0_TU 11_0_TU 11_1_TU 12_0_TU 12_1_TU 13_1_TU 14_1_TU; do
	echo "$subclone gain"
	finalMutFileUnordered=$outDir/${subclone}_gain_TEMP.txt
	finalMutFile=$outDir/${subclone}_gain.vcf

	# sort
	cat $variantListDir/node_${subclone}_gain.txt | sort -k1,1 -k2n > $finalMutFileUnordered

	# add header and mutations to create final vcf file
	echo -e $vcfHeader > $finalMutFile
	cat $finalMutFileUnordered >> $finalMutFile

	echo "make bed file"
	# determine bed file which contains all gain regions that this subclone has
	allGainRegionsBed=$outDir/${subclone}_gain.bed
	cat $cnvListDir/cnv_gains*${subclone}* | sort -k1,1 -k2n | uniq > $outDir/${subclone}_gain_TEMP.bed

	# create final bed file as intersection between the targeted regions and the regions which are gained
	$intersectBed -a $bed_file -b $outDir/${subclone}_gain_TEMP.bed > $allGainRegionsBed
done


