#!/bin/bash

function loc_realign()
{
	if [ -z "$1" -o -z "$2" -o -z "$3" ]; then
		echo "Usage: loc_realign <bamFile> <outputBam> <fn_genome> [ --create-index --check-fexist ]"
	fi
	
	local bamFile=$1; shift
	local outBam=$1; shift
	local fn_genome=$1; shift

	local check_fexist=false
	local create_index=false
	while [ 0 -lt $# ]; do
		if [ $1 == "--create-index" ]; then
			create_index=true
		elif [ $1 == "--check-fexist" ]; then
			check_fexist=true
		fi
		shift
	done
	

        if $check_fexist && [ -f ${outBam} ]; then
                echo "loc_realign: bam file ${outBam} already exists."
                return 0
        fi
	
	local outputIntervals=${outBam%.bam}.intervals
	local tempOutBam=${outBam%.bam}.TEMP.bam


	if $check_fexist && [ ! -f $outputIntervals ]; then
		# Determine (small) suspicious intervals which are likely in need of realignment
	
		# by default read with mapping qualitiy of 255 are filtered out by gatk. However, bowtie assigns mapping quality of 
		# 255 to all reads with only a single alignment. With this option we modify this mapping quality
		rq_opt="-rf ReassignOneMappingQuality -RMQF 255 -RMQT 254"
	
		command="$GATK -T RealignerTargetCreator -R $fn_genome -I $bamFile -o $outputIntervals $rq_opt"
		echo $command
		eval $command
	fi


	if $check_fexist && [ ! -f $tempOutBam ]; then
		# Running the realigner over those intervals 

		# by default read with mapping qualitiy of 255 are filtered out by gatk. However, bowtie assigns mapping quality of 
		# 255 to all reads with only a single alignment. With this option we modify this mapping quality
		rq_opt="-rf ReassignOneMappingQuality -RMQF 255 -RMQT 254"
	
		command="$GATK -T IndelRealigner -R $fn_genome -I $bamFile -targetIntervals $outputIntervals -o $tempOutBam $rq_opt"
		echo $command
		eval $command
	fi


	if $check_fexist && [ ! -f $outBam ]; then
		# Ensure that all mate-pair information is in sync between each read and it's mate pair
		echo "$picard_tools/FixMateInformation.jar  INPUT=$tempOutBam  OUTPUT=$outBam  SORT_ORDER=coordinate  CREATE_INDEX=$create_index "
		$picard_tools/FixMateInformation.jar	INPUT=$tempOutBam \
							OUTPUT=$outBam \
							SORT_ORDER=coordinate \
							CREATE_INDEX=$create_index

		if $create_index; then
			mv ${outBam%.bam}.bai ${outBam}.bai
		fi
	fi
}
