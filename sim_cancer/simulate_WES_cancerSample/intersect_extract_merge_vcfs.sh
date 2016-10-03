#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`
source `find ${gitDir} -name paths.sh`
compute_overlap=`find ${gitDir} -name compute_overlap`

vcfDir=$1

if [ -z "$vcfDir" -o -z "$compute_overlap" ]; then
	echo "No vcfDir or compute_overlap script found"
	exit 0
fi

cd $vcfDir

	normal_all=NO_all.vcf
	normal_smmr=NO_smmr.vcf
		
	tumor_all=TU_all.vcf
	tumor_smmr=TU_smmr.vcf
	
	## generate the intersections
	if false; then 
		
		
		for file1 in $normal_smmr $tumor_all $tumor_smmr
		do
			for file2 in $normal_all
			do
				command="$compute_overlap $file1 $file2 private_${file1%.vcf}_minus_${file2} private_${file2%.vcf}_minus_${file1} overlap_${file1%.vcf}_${file2}"
				echo $command
				eval $command
			done
		done
		
		for file in private_${tumor_all%.vcf}_minus_${normal_all} private_${tumor_smmr%.vcf}_minus_${normal_all}
		do
			command="$compute_overlap $file $normal_smmr private_${file%.vcf}_minus_${normal_smmr} private_${normal_smmr%.vcf}_minus_${file} overlap_${file%.vcf}_${normal_smmr}"
			echo $command
			eval $command
		done
		
		
		file1=private_private_${tumor_all%.vcf}_minus_${normal_all%.vcf}_minus_${normal_smmr}
		file2=private_private_${tumor_smmr%.vcf}_minus_${normal_all%.vcf}_minus_${normal_smmr}
		private_file1=private_private_private_${tumor_all%.vcf}_minus_${normal_all%.vcf}_minus_${normal_smmr%.vcf}_minus_private_private_${tumor_smmr%.vcf}_minus_${normal_all%.vcf}_minus_${normal_smmr}
		private_file2=private_private_private_${tumor_smmr%.vcf}_minus_${normal_all%.vcf}_minus_${normal_smmr%.vcf}_minus_private_private_${tumor_all%.vcf}_minus_${normal_all%.vcf}_minus_${normal_smmr}
		overlap_file1_file2=overlap_private_private_${tumor_all%.vcf}_minus_${normal_all%.vcf}_minus_${normal_smmr%.vcf}_private_private_${tumor_smmr%.vcf}_minus_${normal_all%.vcf}_minus_${normal_smmr}
		
		command="$compute_overlap $file1 $file2 $private_file1 $private_file2 $overlap_file1_file2"
		echo $command
		eval $command
	fi
	
	
	
	
	## extract the percentages and merge the vcf files
	if true; then
		NO_smmr=NO_smmr.vcf
		NO_all_only=private_NO_all_minus_NO_smmr.vcf
		cat $NO_smmr | grep ^# > ${NO_smmr%.vcf}_Header.txt
	
		TU_all_only_only_only=private_private_private_${tumor_all%.vcf}_minus_${normal_all%.vcf}_minus_${normal_smmr%.vcf}_minus_private_private_${tumor_smmr%.vcf}_minus_${normal_all%.vcf}_minus_${normal_smmr}
		TU_smmr_only_only=private_private_${tumor_smmr%.vcf}_minus_${normal_all%.vcf}_minus_${normal_smmr}
		cat $TU_smmr_only_only | grep ^# > ${TU_smmr_only_only%.vcf}_Header.txt
	
		for percentage in 0 0.25 0.33 0.5
		do
			echo "Percentage: $percentage"
			mkdir -p final_vcfs_${percentage}
	
			## The normal
			cat $NO_all_only | grep -v ^# | awk "BEGIN {srand()} { if (rand() <= $percentage) print }" > ./final_vcfs_${percentage}/${NO_all_only%.vcf}_${percentage}_woHeader.txt
	
			echo "cat ./final_vcfs_${percentage}/${NO_all_only%.vcf}_${percentage}_woHeader.txt $NO_smmr | grep -v ^# > ./final_vcfs_${percentage}/NO_${percentage}_unheadered_unsorted.vcf"
			cat ./final_vcfs_${percentage}/${NO_all_only%.vcf}_${percentage}_woHeader.txt $NO_smmr | grep -v ^# > ./final_vcfs_${percentage}/NO_${percentage}_unheadered_unsorted.vcf
	
			echo "cat ./final_vcfs_${percentage}/NO_${percentage}_unheadered_unsorted.vcf | sort -k1,1 -k2n > ./final_vcfs_${percentage}/NO_${percentage}_unheadered.vcf"
			cat ./final_vcfs_${percentage}/NO_${percentage}_unheadered_unsorted.vcf | sort -k1,1 -k2n > ./final_vcfs_${percentage}/NO_${percentage}_unheadered.vcf
	
			echo "cat ${NO_smmr%.vcf}_Header.txt ./final_vcfs_${percentage}/NO_${percentage}_unheadered.vcf > ./final_vcfs_${percentage}/NO_${percentage}_final.vcf"
			cat ${NO_smmr%.vcf}_Header.txt ./final_vcfs_${percentage}/NO_${percentage}_unheadered.vcf > ./final_vcfs_${percentage}/NO_${percentage}_final.vcf
	
			## The tumor
			cat $TU_all_only_only_only | grep -v ^# | awk "BEGIN {srand()} { if (rand() <= $percentage) print }" > ./final_vcfs_${percentage}/${TU_all_only_only_only%.vcf}_${percentage}_woHeader.txt
			
			echo "cat ./final_vcfs_${percentage}/${TU_all_only_only_only%.vcf}_${percentage}_woHeader.txt $TU_smmr_only_only | grep -v ^# > ./final_vcfs_${percentage}/TU_${percentage}_unheadered_unsorted.vcf"
			cat ./final_vcfs_${percentage}/${TU_all_only_only_only%.vcf}_${percentage}_woHeader.txt $TU_smmr_only_only | grep -v ^# > ./final_vcfs_${percentage}/TU_${percentage}_unheadered_unsorted.vcf
	
			echo "cat ./final_vcfs_${percentage}/TU_${percentage}_unheadered_unsorted.vcf | sort -k1,1 -k2n > ./final_vcfs_${percentage}/TU_${percentage}_unheadered.vcf"
			cat ./final_vcfs_${percentage}/TU_${percentage}_unheadered_unsorted.vcf | sort -k1,1 -k2n > ./final_vcfs_${percentage}/TU_${percentage}_unheadered.vcf
	
			echo "cat ${TU_smmr_only_only%.vcf}_Header.txt ./final_vcfs_${percentage}/TU_${percentage}_unheadered.vcf > ./final_vcfs_${percentage}/TU_${percentage}_final.vcf"
			cat ${TU_smmr_only_only%.vcf}_Header.txt ./final_vcfs_${percentage}/TU_${percentage}_unheadered.vcf > ./final_vcfs_${percentage}/TU_${percentage}_final.vcf
		done
	fi
cd -
