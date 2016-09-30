#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`

listDir=$1

if [ -z "$1" ]; then
	echo "Usage: $0 <listDir>"
	exit 0
fi

# callers for venn diagram: 
declare -a callers_venn=('joint' 'deep' 'gatkUG' 'somSniper' 'sinvict') # top 5 callers

for caller_idx in `seq 1 ${#callers_venn[@]}`
do
	wc -l $listDir/${callers_venn[$caller_idx-1]}.txt 
done

already_done="---"
for caller_idx1 in `seq 1 ${#callers_venn[@]}`
do
	for caller_idx2 in `seq 1 ${#callers_venn[@]}`
	do
		if [ "${callers_venn[$caller_idx1-1]}" == "${callers_venn[$caller_idx2-1]}" ]; then
			continue
		fi
		if [[ "$already_done" == *"${callers_venn[$caller_idx2-1]}"* ]]; then
			continue
		fi
		overlap_file=$listDir/overlap_${callers_venn[$caller_idx1-1]}_${callers_venn[$caller_idx2-1]}.txt
		if [ ! -f $overlap_file ]; then
			command="java SNV_Overlapper_Ordered $listDir/${callers_venn[$caller_idx1-1]}.txt $listDir/${callers_venn[$caller_idx2-1]}.txt $overlap_file $listDir/private_${callers_venn[$caller_idx1-1]}_${callers_venn[$caller_idx2-1]}.txt "
			echo $command
			eval $command
		fi
		wc -l $overlap_file 
	done
	already_done="${already_done}_${callers_venn[$caller_idx1-1]}"
done

already_done="---"
for caller_idx1 in `seq 1 ${#callers_venn[@]}`
do
	already_done2="---"
	for caller_idx2 in `seq 1 ${#callers_venn[@]}`
	do
		if [ "${callers_venn[$caller_idx1-1]}" == "${callers_venn[$caller_idx2-1]}" ]; then
			continue
		fi
		if [[ "$already_done2" == *"${callers_venn[$caller_idx2-1]}"* ]]; then
			continue
		fi
		if [[ "$already_done" == *"${callers_venn[$caller_idx1-1]}"* ]]; then
			continue
		fi
		if [[ "$already_done2" == *"${callers_venn[$caller_idx1-1]}"* ]]; then
			continue
		fi
		if [[ "$already_done" == *"${callers_venn[$caller_idx2-1]}"* ]]; then
			continue
		fi
		for caller_idx3 in `seq 1 ${#callers_venn[@]}`
		do

			if [ "${callers_venn[$caller_idx1-1]}" == "${callers_venn[$caller_idx3-1]}" -o "${callers_venn[$caller_idx2-1]}" == "${callers_venn[$caller_idx3-1]}" ]; then
				continue
			fi
			if [[ "$already_done" == *"${callers_venn[$caller_idx3-1]}"* ]]; then
				continue
			fi
			if [[ "$already_done2" == *"${callers_venn[$caller_idx3-1]}"* ]]; then
				continue
			fi
			overlap_file=$listDir/overlap_${callers_venn[$caller_idx1-1]}_overlap_${callers_venn[$caller_idx2-1]}_${callers_venn[$caller_idx3-1]}.txt
			if [ ! -f $overlap_file ]; then
				command="java SNV_Overlapper_Ordered $listDir/${callers_venn[$caller_idx1-1]}.txt $listDir/overlap_${callers_venn[$caller_idx2-1]}_${callers_venn[$caller_idx3-1]}.txt $overlap_file $listDir/private_${callers_venn[$caller_idx1-1]}_overlap_${callers_venn[$caller_idx2-1]}_${callers_venn[$caller_idx3-1]}.txt "
				echo $command
				eval $command
			fi
			wc -l $overlap_file
		done
		already_done2="${already_done2}_${callers_venn[$caller_idx2-1]}"
	done
	already_done="${already_done}_${callers_venn[$caller_idx1-1]}"
done


for file in $listDir/overlap_${callers_venn[0]}_overlap_${callers_venn[1]}_overlap_${callers_venn[2]}_${callers_venn[3]}.txt $listDir/overlap_${callers_venn[0]}_overlap_${callers_venn[1]}_overlap_${callers_venn[2]}_${callers_venn[4]}.txt $listDir/overlap_${callers_venn[0]}_overlap_${callers_venn[1]}_overlap_${callers_venn[3]}_${callers_venn[4]}.txt $listDir/overlap_${callers_venn[0]}_overlap_${callers_venn[2]}_overlap_${callers_venn[3]}_${callers_venn[4]}.txt $listDir/overlap_${callers_venn[1]}_overlap_${callers_venn[2]}_overlap_${callers_venn[3]}_${callers_venn[4]}.txt
do
	if [ ! -f $file ]; then
		caller1=`basename $file | tr '_' ' ' | awk '{print $2}'`
		caller2=`basename $file | tr '_' ' ' | awk '{print $4}'`
		caller3=`basename $file | tr '_' ' ' | awk '{print $6}'`
		caller4=`basename ${file%.txt} | tr '_' ' ' | awk '{print $7}'`
		command="java SNV_Overlapper_Ordered $listDir/$caller1.txt $listDir/overlap_${caller2}_overlap_${caller3}_${caller4}.txt $file $listDir/private_${caller1}_overlap_${caller2}_overlap_${caller3}_${caller4}.txt"
		echo $command
		eval $command
	fi
	wc -l $file
done

overlap_file=$listDir/overlap_${callers_venn[0]}_overlap_${callers_venn[1]}_overlap_${callers_venn[2]}_overlap_${callers_venn[3]}_${callers_venn[4]}.txt
if [ ! -f $overlap_file ]; then
	command="java SNV_Overlapper_Ordered $listDir/${callers_venn[0]}.txt $listDir/overlap_${callers_venn[1]}_overlap_${callers_venn[2]}_overlap_${callers_venn[3]}_${callers_venn[4]}.txt $overlap_file $listDir/private_${callers_venn[0]}_overlap_${callers_venn[1]}_overlap_${callers_venn[2]}_overlap_${callers_venn[3]}_${callers_venn[4]}.txt"
	echo $command
	eval $command
fi
wc -l $overlap_file 

