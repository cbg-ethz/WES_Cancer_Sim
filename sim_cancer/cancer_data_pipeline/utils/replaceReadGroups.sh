#!/bin/bash

function replaceReadGroups()
{

	if [ -z "$1" -o -z "$2" -o -z "$3" ]; then
		echo "Usage: replaceReadGroups <bamFile> <outputBam> <RGID> <RGLB> <RGPL> <RGPU> <RGSM> [ --index-bam <false/true> ]"
	fi
	
	local bamFile=$1; shift
	local outBam=$1; shift
	local RGID=$1; shift
	local RGLB=$1; shift
	local RGPL=$1; shift
	local RGPU=$1; shift
	local RGSM=$1; shift
	local check_fexist=false

	while [ 0 -lt $# ]; do
		if [ $1 == "--check-fexist" ]; then
			check_fexist=true
		else
			echo replaceReadGroups: did not understand arg $1
			exit 0
		fi;
		shift
	done
	
	
	if [  $check_fexist -a -f ${outBam} ]; then
		echo "replaceReadGroups: bam file ${outBam} already exists."
		return 0
	fi
		
	echo "$picard_tools/AddOrReplaceReadGroups.jar INPUT=$bamFile OUTPUT=$outBam RGID=$RGID RGLB=$RGLB RGPL=$RGPL RGPU=$RGPU RGSM=$RGSM"
	$picard_tools/AddOrReplaceReadGroups.jar INPUT=$bamFile OUTPUT=$outBam RGID=$RGID RGLB=$RGLB RGPL=$RGPL RGPU=$RGPU RGSM=$RGSM

}
