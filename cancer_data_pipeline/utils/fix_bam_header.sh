#!/bin/bash

function fix_bam_header()
{

	if [ -z "$1" -o -z "$2" -o -z "$3" ]; then
		echo "Usage: $0 <bamFiles> <extension> <outputDir> [ --run-opt <submit> ]"
		exit -1
	fi
	
	local bamFiles=$1
	local ext=$2
	local outdir=$3
	
	for i in {1..3}
	do
		shift;
	done
	
	local run_opt="local"
	
	while [ 1 -lt $# ]; do
		if [ $1 == "--run-opt" ]; then
			shift;
			run_opt=$1;
		else
			echo did not understand arg $1
			shift
		fi;
	done
	
	
	mkdir -p $outdir
	

	
	for fn_bam in `echo $bamFiles`; do 
		
		local fn_result=$outdir/`basename ${fn_bam%$ext}`.bam
	
		if [ -f $fn_result ]; then 
			echo "bam file $fn_result already exists."
			continue; 
		elif [ -h $fn_result ]; then # symbolic link exists
			echo "bam file $fn_result already exists."
                        continue;
		else
			echo "bam file $fn_result does not yet exist."			
		fi
				

		if [ $run_opt == "submit" ]; then
			echo "submit $fn_bam"
			local log_dir=$outdir/logs/
			mkdir -p $log_dir
			local fn_out=${log_dir}reheader_`basename ${fn_bam}`.o
			local fn_err=${log_dir}reheader_`basename ${fn_bam}`.e
			echo "" > $fn_err
			echo "" > $fn_out
			num_jobs=1000
			while [ "$num_jobs" -gt "15" ]
			do
				sleep 30
				num_jobs=`qstat | grep reheader | wc -l`
				echo num_jobs: $num_jobs
			done
			queue="-q regular.q@bs-dsvr02,regular.q@bs-dsvr04,regular.q@bs-dsvr05,regular.q@bs-dsvr08"
			echo "echo ${fn_bam} > $fn_out; hostname >> $fn_out; source `find ${gitDir} -name fix_bam_header.sh` ; source `find ${gitDir} -name paths.sh`; source `find ${gitDir} -name utils.sh`;  fix_bam_header $fn_bam $ext $outdir " | qsub $queue -o $fn_out -e $fn_err -N reheader_`basename ${fn_bam}` -cwd -V
		else
	

			unique_header=$outdir/`basename ${fn_bam}`_unique_header.txt
			if [ ! -f $unique_header ]; then
				echo "$samtools view -H $fn_bam | uniq > $unique_header"
				$samtools view -H $fn_bam | uniq > $unique_header
			else
				echo "$unique_header already exists."
			fi
	
		
			if [ ! -f $fn_result ]; then
				echo "$picard_tools/ReplaceSamHeader.jar INPUT=$fn_bam HEADER=$unique_header OUTPUT=$fn_result"	
			        $picard_tools/ReplaceSamHeader.jar INPUT=$fn_bam HEADER=$unique_header OUTPUT=$fn_result
			else
				echo "$fn_result already exists."
			fi
	

			if [ ! -f ${fn_result}.bai ]; then
			        echo "$samtools index $fn_result"
			        samtools index $fn_result
			else
				echo "${fn_result}.bai already exists."
			fi

		fi

	done
}
