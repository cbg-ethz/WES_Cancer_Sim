#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`
source `find ${gitDir} -name paths.sh`
tools_dir=${gitDir}
source `find ${tools_dir} -name genome.sh`

vcf_dir=$1 # this folder has to contain the following files: 
## 0_NO_final.vcf 1_NO_final.vcf
## 
## 7_0_TU_final.vcf 7_1_TU_final.vcf   8_0_TU_final.vcf 8_1_TU_final.vcf    9_0_TU_final.vcf 9_1_TU_final.vcf   10_0_TU_final.vcf 10_1_TU_final.vcf                                                                  11_0_TU_final.vcf 11_1_TU_final.vcf   12_0_TU_final.vcf 12_1_TU_final.vcf   13_0_TU_final.vcf 13_1_TU_final.vcf   14_0_TU_final.vcf 14_1_TU_final.vcf   

out_dir=$2

if [ -z "$vcf_dir" -o -z "$out_dir" ]; then
	echo "vcf_dir or out_dir not provided"
	exit 0
fi


mkdir -p $out_dir

# create a vcf file without header
for vcf in $vcf_dir/[0-9]*final.vcf; do
	outfile=$out_dir/`basename ${vcf%.vcf}`_woHeader.txt
	if [ ! -f $outfile ]; then
		cat $vcf | grep -v ^# > $outfile
	fi
done

# create intersections to get vcf files from intermediate subclones
for diploidCopy in 0 1; do
	for cloneCnt in 7 9 11 13; do
		neighbor=$((cloneCnt+1))

		file1=$out_dir/${cloneCnt}_${diploidCopy}_TU_final_woHeader.txt
		file2=$out_dir/${neighbor}_${diploidCopy}_TU_final_woHeader.txt
		overlapFile1File2=$out_dir/overlap_`basename ${file1%_final_woHeader.txt}`_`basename ${file2%_final_woHeader.txt}`.txt
		privateFile1_minusFile2=$out_dir/private_`basename ${file1%_final_woHeader.txt}`_minus_`basename ${file2%_final_woHeader.txt}`.txt
		if [ ! -f $overlapFile1File2 ]; then
			command="java SNV_Overlapper_Ordered_PlusAlleles $file1 $file2 $overlapFile1File2 $privateFile1_minusFile2"
			echo $command
			eval $command
		fi

		file1=$out_dir/${neighbor}_${diploidCopy}_TU_final_woHeader.txt
		file2=$out_dir/${cloneCnt}_${diploidCopy}_TU_final_woHeader.txt
		overlapFile1File2=$out_dir/overlap_`basename ${file1%_final_woHeader.txt}`_`basename ${file2%_final_woHeader.txt}`.txt
		privateFile1_minusFile2=$out_dir/private_`basename ${file1%_final_woHeader.txt}`_minus_`basename ${file2%_final_woHeader.txt}`.txt
		if [ ! -f $overlapFile1File2 ]; then
			command="java SNV_Overlapper_Ordered_PlusAlleles $file1 $file2 $overlapFile1File2 $privateFile1_minusFile2"
			echo $command
			eval $command
		fi
	done

	for cloneCnt in 7 11; do
		directNeighbor=$((cloneCnt+1))
		oneFurther=$((cloneCnt+2))
		twoFurther=$((cloneCnt+3))

		file1=$out_dir/overlap_${cloneCnt}_${diploidCopy}_TU_${directNeighbor}_${diploidCopy}_TU.txt
		file2=$out_dir/overlap_${oneFurther}_${diploidCopy}_TU_${twoFurther}_${diploidCopy}_TU.txt
		overlapFile1File2=$out_dir/overlap_${cloneCnt}_${diploidCopy}_TU_${directNeighbor}_${diploidCopy}_TU_${oneFurther}_${diploidCopy}_TU_${twoFurther}_${diploidCopy}_TU.txt
		privateFile1_minusFile2=$out_dir/private_${cloneCnt}_${diploidCopy}_TU_${directNeighbor}_${diploidCopy}_TU_${oneFurther}_${diploidCopy}_TU_${twoFurther}_${diploidCopy}_TU.txt
		if [ ! -f $overlapFile1File2 ]; then
                        command="java SNV_Overlapper_Ordered_PlusAlleles $file1 $file2 $overlapFile1File2 $privateFile1_minusFile2"
                        echo $command
                        eval $command
                fi

		file1=$out_dir/overlap_${oneFurther}_${diploidCopy}_TU_${twoFurther}_${diploidCopy}_TU.txt
		file2=$out_dir/overlap_${cloneCnt}_${diploidCopy}_TU_${directNeighbor}_${diploidCopy}_TU.txt
		overlapFile1File2=$out_dir/overlap_${oneFurther}_${diploidCopy}_TU_${twoFurther}_${diploidCopy}_TU_${cloneCnt}_${diploidCopy}_TU_${directNeighbor}_${diploidCopy}_TU.txt
		privateFile1_minusFile2=$out_dir/private_${oneFurther}_${diploidCopy}_TU_${twoFurther}_${diploidCopy}_TU_${cloneCnt}_${diploidCopy}_TU_${directNeighbor}_${diploidCopy}_TU.txt
		if [ ! -f $overlapFile1File2 ]; then
                        command="java SNV_Overlapper_Ordered_PlusAlleles $file1 $file2 $overlapFile1File2 $privateFile1_minusFile2"
                        echo $command
                        eval $command
                fi	
	done

	file1=$out_dir/overlap_7_${diploidCopy}_TU_8_${diploidCopy}_TU_9_${diploidCopy}_TU_10_${diploidCopy}_TU.txt
	file2=$out_dir/overlap_11_${diploidCopy}_TU_12_${diploidCopy}_TU_13_${diploidCopy}_TU_14_${diploidCopy}_TU.txt
	overlapFile1File2=$out_dir/overlap_7891011121314_${diploidCopy}_TU.txt
	privateFile1_minusFile2=$out_dir//private_7891011121314_${diploidCopy}_TU.txt
	if [ ! -f $overlapFile1File2 ]; then
		command="java SNV_Overlapper_Ordered_PlusAlleles $file1 $file2 $overlapFile1File2 $privateFile1_minusFile2"
		echo $command
		eval $command
	fi

	file1=$out_dir/overlap_11_${diploidCopy}_TU_12_${diploidCopy}_TU_13_${diploidCopy}_TU_14_${diploidCopy}_TU.txt
	file2=$out_dir/overlap_7_${diploidCopy}_TU_8_${diploidCopy}_TU_9_${diploidCopy}_TU_10_${diploidCopy}_TU.txt
	overlapFile1File2=$out_dir/overlap_1112131478910_${diploidCopy}_TU.txt
	privateFile1_minusFile2=$out_dir/private_1112131478910_${diploidCopy}_TU.txt
	if [ ! -f $overlapFile1File2 ]; then
		command="java SNV_Overlapper_Ordered_PlusAlleles $file1 $file2 $overlapFile1File2 $privateFile1_minusFile2"
		echo $command
		eval $command
	fi
done

exit 0


# for fn_vcf in $out_dir/[0-9]*.vcf
# do
# 	fn_ref_out=${fn_vcf%.vcf}.fa
# 	fn_map=${fn_vcf%.vcf}_mapping.csv
# 	command="$create_ref_from_vcf $ref $fn_vcf $fn_ref_out $fn_map >> log.log 2>&1 " 
# 	echo $command
# 	eval $command
# done

