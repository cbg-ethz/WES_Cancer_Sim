#!/bin/bash

function sort_bam()
{
	source `find ${gitDir} -name paths.sh`

	if [ -z "$1" -o -z "$2" -o -z "$3" ]; then
		echo "Usage: sort_bam <bamFiles> <outputDir> <statusDir> [ --run-opt <submit> --bam-ext <ext> --human-karyotypic <false/true> ]"
	fi
	
	local bamFileS=$1; shift; 
	local outdir=$1; shift; 
	local statusDir=$1; shift;
	local ext=.bam	
	local human_karyotypic=false

	while [ 0 -lt $# ]; do
		if [ $1 == "--bam-ext" ]; then
			shift;
			ext=$1;
		elif [ $1 == "--human-karyotypic" ]; then
			human_karyotypic=true
		else
			echo did not understand arg $1
		fi;
		shift
	done
	
	echo "bam file extension is $ext"
	echo "sort bam file in karyotypic order (only for human data): $human_karyotypic"
	
	mkdir -p $outdir
		
	for bamFile in `echo $bamFileS `; do 
		
		local fn_result=$outdir/`basename ${bamFile%.bam}`.sorted.bam
	
		if $human_karyotypic; then
			if [ -s ${bamFile%.bam}.karyotypic.bam ]; then
				echo "bam file ${bamFile%.bam}.karyotypic.bam already exists."	
				continue;
			fi
		else
			if [ -s ${fn_result} ]; then
				echo "bam file ${fn_result} already exists."
				continue;
			fi
		fi
	
		if $human_karyotypic; then
			echo "$picard_tools/ReorderSam.jar INPUT=$bamFile OUTPUT=${bamFile%.bam}.karyotypic.bam REFERENCE=${humanRefKaryo}"
			$picard_tools/ReorderSam.jar INPUT=$bamFile OUTPUT=${bamFile%.bam}.karyotypic.bam REFERENCE=${humanRefKaryo}
		else
			echo "$picard_tools/SortSam.jar INPUT=$bamFile OUTPUT=$fn_result SORT_ORDER=coordinate"
			$picard_tools/SortSam.jar INPUT=$bamFile OUTPUT=$fn_result SORT_ORDER=coordinate
		fi
	done
}
