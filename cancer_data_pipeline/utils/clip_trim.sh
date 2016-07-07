#!/bin/bash

function clip_trim()
{

	if [ -z $1 ]; then
		echo "Usage: $0 <input_dir> <output_dir> <clip> <trim> <AdapterFilePrefix> <nr_threads> [ --no-single <false/true> --run-opt <submit/check> --tag1 <_R1> --tag2 <_R2> --score <12> ]"
	fi
	
	local openbis_dir=$1
	local output_dir=$2
	local clip=$3
	local trim=$4
	local AdapterFilePrefix=$5
	local p=$6 #threads
	
	for i in {1..6}
	do
		shift;  
	done
	
	local no_single=false
	local run_opt=""
	local tag1=_R1
	local tag2=_R2
	local score=7
	
	local i=1
	while [ $i -lt $# ]; do 
		if [ $1 == "--no-single" ]; then 
			shift; 
			no_single=$1; 
		elif [ $1 == "--run-opt" ]; then 
			shift; 
			run_opt=$1; 
		elif [ $1 == "--tag1" ]; then 
			shift; 
			tag1=$1; 
		elif [ $1 == "--tag2" ]; then 
			shift; 
			tag2=$1; 
		elif [ $1 == "--score" ]; then
			shift;
			score=$1
		else
			echo did not understand arg $1
			shift;
		fi; 
	done
	
	echo "Using the following options:"
	echo "clip single read fastq files = $no_single"
	if [ -z "$run_opt" -o "$run_opt" == "local" ]; then
		echo "Running the job locally"
	else
		echo "Submitting the job to the grid"
	fi
	echo "Tags for the forward/reverse reads are $tag1 and $tag2"
	echo "The score parameter for clipping is set to $score"
	
	
	mkdir -p $output_dir
	
	
	################################################################################
	# if arg7 not specified, then submit all jobs to cluster
	################################################################################
	if [ "$run_opt" == "submit" ]; then
		echo submit $run_opt
	
		#####################
		# stats
		local cnt_todo=0
		local cnt_done=0
		for f in `find $openbis_dir -name \*.tsv`; do 				       
			command="$0 `dirname $f` $output_dir $clip $trim $AdapterFilePrefix $p --no-single $no_single --run-opt check --tag1 $tag1 --tag2 $tag2 --score $score"
			if [ -z "`$command | grep TODO`" ]; then 
				cnt_done=$((cnt_done+1)) 
			else
				cnt_todo=$((cnt_todo+1))
			fi
		done
		echo clipping/trimming done for $cnt_done tsv-files todo: $cnt_todo
	
		if [ "$cnt_todo" == "0" ]; then 
			exit 0
		fi
		
		#####################
		# run it 
		local cnt=0
		local log_dir=$output_dir/logs/
		local dir_=`pwd`
		mkdir -p $log_dir
		for f in `find $openbis_dir -name \*.tsv`; do
			cnt=$((cnt + 1))
			local fn_out=${log_dir}clip_`basename $f`.o 
			local fn_err=${log_dir}clip_`basename $f`.e
			echo "" > $fn_err
			echo "" > $fn_out
			echo "clip number (cnt) = $cnt" >> $fn_out
	
			local command="$0 `dirname $f` $output_dir $clip $trim $AdapterFilePrefix $p  --no-single $no_single --run-opt check --tag1 $tag1 --tag2 $tag2 --score $score"
			if [ -z "`$command | grep TODO`" ]; then 
				continue
			fi
			echo running for file $f
			local num_jobs=1000
			while [ "$num_jobs" -gt "15" ]
			do
				sleep 10
				num_jobs=`qstat | grep clip | wc -l`
				echo num_jobs: $num_jobs
			done
			local queue="-q regular.q"
			#queue="-q mpi04-ht.q" # bewi												 
			echo "cd $dir_ ; echo $f > $fn_out; hostname >> $fn_out; $0 `dirname $f` $output_dir $clip $trim $AdapterFilePrefix $p --no-single $no_single --tag1 $tag1 --tag2 $tag2 --score $score "| qsub $queue -o $fn_out -e $fn_err -N clip$cnt -cwd -V -pe make $p
		done
		exit 0
	fi
	################################################################################
	
	
	
	for tsv in `find $openbis_dir -name \*.tsv`; do
		
		local sample=`extract $tsv EXTERNAL_SAMPLE_NAME`
		local protocol=`extract $tsv END_TYPE`
	
		local dd=`dirname $tsv`
		local outdir=$output_dir/fastqs/`basename ${tsv%.tsv}`/
		mkdir -p $outdir
		mkdir -p $outdir/paired/
	
		if [ ! "$run_opt" == "check" ]; then 
			ln -s $tsv $outdir/paired/`basename $tsv` 
		fi
	
		local fn_stats=$output_dir/stats/`basename ${tsv%.tsv}`.stats
		mkdir -p `dirname $fn_stats`
	
		for f in `ls $dd | grep $tag1`; do
			local fastq1=$dd/$f
			if [ -z $fastq1 ]; then
				echo did not find fastq file in $dd
				exit -1
			fi
			local fastq2=`echo $fastq1 | sed "s/$tag1/$tag2/"`
	
			local trim_report=$outdir/logs/`basename ${tsv%.tsv}`_cliptrim.log
			mkdir -p `dirname $trim_report`
	
			local opts1="-phred33 -trimlog $trim_report "
			local opts=" "
			if [ $trim == 1 ];
			then
				local minqual=2
				opts="SLIDINGWINDOW:4:$minqual LEADING:$minqual TRAILING:$minqual MINLEN:50 "
			else
				opts="$opts MINLEN:50 "
			fi
	
			if [ $protocol == "PAIRED_END" ]; then
	
				if [ ! -f ${AdapterFilePrefix}_PE.fa ]; then
					echo "Adapter file ${AdapterFilePrefix}_PE.fa does not exist."
					exit -1
				fi
				count_fastq $fn_stats ${fastq1} ${fastq2}
	
				if [ $clip == 1 ];
				then
					#ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
					seedmismatches=1
					palindrom=30 # recommended by trimmomatic ## this is only for palindrome mode (i.e. when adpater file is in special format)
					##score=7 # perfect match of 12 bases
					##score=15 # perfect match of 25 bases
					min_adapt_len=1 # recommended by trimmomatic ## this is only for palindrome mode
					keepboth=true # do not discard one read if pair is overlapping ## this is only for palindrome mode 
					opts=" ILLUMINACLIP:${AdapterFilePrefix}_PE.fa:$seedmismatches:$palindrom:$score:$min_adapt_len:$keepboth $opts "
				fi
	
				local PairedOutput1=$outdir/paired/`basename ${fastq1%.fastq.gz}`_cliptrim.fastq.gz
				local PairedOutput2=$outdir/paired/`basename ${fastq2%.fastq.gz}`_cliptrim.fastq.gz
				local UnPairedOutput1=$outdir/unpaired1/`basename ${fastq1%.fastq.gz}`_UnPcliptrim.fastq.gz
				local UnPairedOutput2=$outdir/unpaired2/`basename ${fastq2%.fastq.gz}`_UnPcliptrim.fastq.gz
				mkdir -p $outdir/unpaired1/
				mkdir -p $outdir/unpaired2/
				cat $tsv | sed 's/END_TYPE\tPAIRED_END/END_TYPE\tSINGLE_READ/' > $outdir/unpaired1/`basename ${tsv%.tsv}`_UnP1.tsv
				cat $tsv | sed 's/END_TYPE\tPAIRED_END/END_TYPE\tSINGLE_READ/' > $outdir/unpaired2/`basename ${tsv%.tsv}`_UnP2.tsv
	
				#if [ ! -f $trim_report -a ! -s $PairedOutput1 -a ! -s $PairedOutput2 -a ! -f $UnPairedOutput1 -a ! -f $UnPairedOutput2 ]; then
				if [ ! -f asdfasdf ]; then
					if [ "$run_opt" == "check" ] 
					then
						echo TODO
					else
						echo "The clip_trim option is set to $opts"
						$trimmomatic PE	 -threads $p \
								$opts1 \
								$fastq1 \
								$fastq2 \
								$PairedOutput1 \
								$UnPairedOutput1 \
								$PairedOutput2 \
								$UnPairedOutput2 \
								$opts
	
						#######################################
						# check and stats
						echo "after clip/trim" >> $fn_stats
						count_fastq $fn_stats $PairedOutput1 $PairedOutput2 $UnPairedOutput1 $UnPairedOutput2
						check_gzip $fn_stats $PairedOutput1 $PairedOutput2 $UnPairedOutput1 $UnPairedOutput2
						#######################################
					fi
				else
					echo "Clipping and/or trimming was already done for `basename $fastq1` and `basename $fastq2`."
					continue
				fi
	
			else ## single end fastq files
	
				if $no_single ; then
					continue;
				fi
	
				if [ ! -f ${AdapterFilePrefix}_SE.fa ]; then
					echo "Adapter file ${AdapterFilePrefix}_SE.fa does not exist."
					exit -1
				fi
	
				count_fastq $fn_stats ${fastq1}
	
				if [ $clip == 1 ];
				then
					#ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
	                                seedmismatches=1
	                                palindrom=30 # recommended by trimmomatic ## this is only for palindrome mode (i.e. when adpater file is in special format)
	                                ##score=7 # perfect match of 12 bases
	                                ##score=15 # perfect match of 25 bases
	
					opts=" ILLUMINACLIP:${AdapterFilePrefix}_SE.fa:$seedmismatches:$palindrom:$score $opts "
				fi
	
				Output=$outdir/`basename ${fastq1%.fastq.gz}`_cliptrim.fastq.gz
	
				if [ ! -s $trim_report -a ! -s $Output ]; then
					if [ "$run_opt" == "check" ] 
					then
						echo TODO
					else
						echo "The clip_trim option is set to $opts"
						$trimmomatic SE	 -threads $p \
								$opts1 \
								$fastq1 \
								$Output \
								$opts
						#######################################
						# check and stats
						count_fastq $fn_stats $Output
						check_gzip $fn_stats $Output
						#######################################
					fi
				else
					echo "Clipping and/or trimming was already done for `basename $fastq1`."
					continue
				fi
	
	
			fi
		done
	done
}
