#!/bin/bash

function merge_bams()
{

	if [ -z "$1" -o -z "$2" ]; then
		echo "Usage: merge_bams <bamFilesToMerge> <outputBam>"
	fi
	
	local bamFileS=$1
	local outBam=$2
	
	if [ -f $outBam ]; then
		echo "bam file $outBam already exists."	
		exit -1
	else
		echo "$samtools merge $outBam `echo $bamFileS`"
		$samtools merge $outBam `echo $bamFileS`
	fi
		
}
