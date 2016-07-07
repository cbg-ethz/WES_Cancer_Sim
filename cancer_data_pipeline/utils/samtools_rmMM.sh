#!/bin/bash

function samtools_rmMM()
{

	if [ -z $1 ]; then
		echo "Usage: samtools_rmMM <bamFileS>"
	fi
	
	bamFileS=$1
	
	for bam_file in `echo $bamFileS`; do 
		
		bam_result=${bam_file%.bam}.woMM.bam
	
		if [ -f bam_result ]; then 
			echo "bam file $bam_result already exists."
			continue; 
		fi
	
		echo "$samtools view -bh -F 256 $bam_file > $bam_result"
		$samtools view -bh -F 256 $bam_file > $bam_result
	
	done
}
