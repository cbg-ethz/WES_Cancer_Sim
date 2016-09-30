#!/bin/bash

function subsetSam()
{

	if [ -z "$1" -o -z "$2" -o -z "$3" ]; then
		echo "Usage: subsetSam <bamFile> <outputBam> <fraction> [ --check-fexist --my-seed <INT> ]"
	fi
	
	local bamFile=$1; shift
	local outBam=$1; shift
	local fraction=$1; shift
	local check_fexist=false
	local my_seed=null

	while [ 0 -lt $# ]; do
		if [ $1 == "--check-fexist" ]; then
			check_fexist=true
		elif [ $1 == "--my-seed" ]; then
			shift;
			my_seed=$1
		else
			echo did not understand arg $1
		fi;
		shift
	done
	
	if $check_fexist && [ -f ${outBam} ]; then
		echo "subsetSam: bam file ${outBam} already exists."
		return 0
	fi
	
	command="$picard_tools/DownsampleSam.jar INPUT=$bamFile OUTPUT=$outBam RANDOM_SEED=$my_seed PROBABILITY=$fraction"
	echo $command
	eval $command
	
}
