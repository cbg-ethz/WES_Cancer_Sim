#!/bin/bash
SGE_stuff="export SGE_ROOT=/usr/local/grid/soge; export BSSE_SGE=1; export SGE_CELL=bsse; export SGE_ARCH=lx-amd64; export SGE_CLUSTER_NAME=bsse-SoGE; export PATH=/usr/local/grid/soge/bin/lx-amd64/:\$PATH"
submit_host=bs-submit01.ethz.ch

function sub_qstat()
{
	qstat 2>/dev/null >> /dev/null
	local ret=$?
	local is_subhost=false
	if [ "$ret" -eq "0" ]; then 
		is_subhost=true
	fi
	
	if $is_subhost; then 
		echo run qstat locally
		qstat $*
	else
		ssh $submit_host "$SGE_stuff; qstat $* " 2> /dev/null 
	fi

}

function wait_for_jobs()
{
	max_num=0
	tag=$1; shift
	while [ 0 -lt $# ]; do  
		if [ $1 == "--max-num" ]; then
			shift
			max_num=$1
		fi
		shift
	done

	sleeptime=1
	num_jobs=1000
	while [ "$num_jobs" -gt "$max_num" ]
	do
		num_jobs=`sub_qstat | grep ${tag:0:9} | wc -l`

		if [ $num_jobs -eq 0 ]; then 
			break; 
		else
			echo "num_jobs(\"$tag\"): $num_jobs sleep:$sleeptime"
			sleep $sleeptime
			sleeptime=$((sleeptime + 1))
		fi
	done
}

function submit()
{
	#defaults
	################################################################################
	#local queue="mpi04-ht.q@bs-dsvr56.ethz.ch"
	local queue="mpi04-ht.q"
	local work_dir=`pwd`
	local job_tag=SUB_`date +%H:%M:%S_%h_%d_%y`_${RANDOM}
	local log_dir=$HOME/tmp/qsub_logs
	local max_num_jobs=100
	local wait_for_jobs=false
	local wait_tag=SUB
	local fn_out=
	local fn_err=
	
	# parse args
	################################################################################
	local command="$1"; shift; 
	
	while [ 0 -lt $# ]; do  
		if [ $1 == "--queue" ]; then
			shift
			queue=$1
		elif [ $1 == "--tag" ]; then
			shift
			job_tag=$1;
		elif [ $1 == "--log-dir" ]; then
			shift
			log_dir=$1;
		elif [ $1 == "--fn-out" ]; then
			shift
			fn_out=$1
			log_dir=`dirname $fn_out`
		elif [ $1 == "--work-dir" ]; then
			shift
			work_dir=$1;
		elif [ $1 == "--wait-tag" ]; then
			shift
			wait_tag=$1;
			wait_for_jobs=true
		elif [ $1 == "--max-jobs" ]; then
			shift
			max_num_jobs=$1
		fi
		shift
	done


	# sleep before submitting a job if queue is full
	################################################################################
	if $wait_for_jobs; then
		num_jobs=1000
		sleep_time=1
		while [ "$num_jobs" -gt $max_num_jobs ]
		do
			num_jobs=`sub_qstat | grep $wait_tag | wc -l`
			echo "num_jobs: $num_jobs > $max_num_jobs, sleep $sleep_time" 
			sleep $sleep_time
			sleep_time=$((sleep_time + 1))
		done
	fi
	
	# determine if this is a submit host
	################################################################################
	qstat 2>/dev/null >> /dev/null
	local ret=$?
	local is_subhost=false
	
	if [ "$ret" -eq "0" ]; then 
		is_subhost=true
	fi
	if [ -z "$fn_out" ]; then 
		fn_out=$log_dir/${job_tag}_`date +%H.%M.%S_%h_%d_%y`_${RANDOM}.o
	fi
	if [ -z "$fn_err" ]; then 
		fn_err=${fn_out%.o}.e
	fi

	local run_on_subhost="mkdir -p $log_dir; echo > $fn_out; echo > $fn_err "


	
	if $is_subhost; then 
		echo "$run_on_subhost" | /bin/bash
		echo "echo \"cd $work_dir; echo PATH:; echo \$PATH;  hostname; echo $command ;  $command \" | /bin/bash" | qsub -q $queue -o $fn_out -e $fn_err -N $job_tag
	else
		ssh $submit_host "$SGE_stuff; $run_on_subhost; echo \"echo \\\"cd $work_dir; echo PATH:; echo \\\$PATH;  hostname;  $command \\\" | /bin/bash \"  | qsub -q $queue -o $fn_out -e $fn_err -N $job_tag " 2>/dev/null
	fi

}
