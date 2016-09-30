#!/bin/bash

function picard_rmdup()
{
	local outdir=$1
	local ext=$2
	local run_opt=$3
	
	if [ "$run_opt" == "submit" ]; then 
		echo submit $0 
		cnt=0
		log_dir=$outdir/logs/
		dir_=`pwd`
		mkdir -p $log_dir
		for fn_bam in $outdir/*$ext; do
			cnt=$((cnt + 1))
			fn_out=${log_dir}rmdup_job$cnt.o
			fn_err=${log_dir}rmdup_job$cnt.e
			echo "" > $fn_err
			echo "" > $fn_out
			num_jobs=1000
			while [ "$num_jobs" -gt "15" ]
			do
				sleep 10
				num_jobs=`qstat | grep rmdup | wc -l`
				echo num_jobs: $num_jobs
			done
			#queue="-q regular.q" 
			queue="-q mpi04-ht.q" # bewi 
			#echo "cd $dir_ ; echo $fn_bam > $fn_out; hostname >> $fn_out; echo shell:$SHELL >> $fn_out ; $0 $outdir `basename $fn_bam` | qsub $queue -o $fn_out -e $fn_err -N rmdup$cnt"
			echo "cd $dir_ ; echo $fn_bam > $fn_out; hostname >> $fn_out; echo shell:$SHELL >> $fn_out ; $0 $outdir `basename $fn_bam`" | qsub $queue -o $fn_out -e $fn_err -N rmdup$cnt
		done
		exit 0
	fi
	
	for fn_bam in $outdir/*$ext; do
	
		local fn_picard=${fn_bam%bam}picard.bam
		local fn_rmdup=${fn_bam%bam}rmdup.bam
		local fn_log=${fn_bam%bam}rmdup.log
	
		opts="INPUT=$fn_bam"
		opts="$opts OUTPUT=$fn_picard"
		opts="$opts SORT_ORDER=coordinate"
	
		if [ ! -s $fn_picard ]; then 
			echo $picard_tools/FixMateInformation.jar $opts 
			$picard_tools/FixMateInformation.jar $opts 
			ret=$?
			if [ $ret != 0 ]; then 
				echo $picard_tools/FixMateInformation.jar $opts returned $ret 
				exit 0
			fi
			echo
		fi
	
		opts="INPUT=$fn_picard"
		opts="$opts OUTPUT=$fn_rmdup"
		opts="$opts METRICS_FILE=$fn_log"
		opts="$opts REMOVE_DUPLICATES=true"
		opts="$opts CREATE_INDEX=false"
		#opts="$opts QUIET=true"
	
		if [ ! -s $fn_rmdup ]; then 
			echo $picard_tools/MarkDuplicates.jar $opts
			$picard_tools/MarkDuplicates.jar $opts
			ret=$?
			if [ $ret != 0 ]; then 
				echo $picard_tools/MarkDuplicates.jar $opts returned $ret 
				exit 0
			fi
	
			echo
		fi
	done
}
