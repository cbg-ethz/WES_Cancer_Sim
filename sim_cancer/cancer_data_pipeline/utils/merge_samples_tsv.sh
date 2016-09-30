#!/bin/bash

function get_tsv_for_bam()
{
	local fn_bam=$1
	local open_bis=$2
	local ext=$3
	
	local base=`basename $fn_bam`
	base=`echo $base | sed 's/_UnP[12]//'`

	local tsv=`find $open_bis -name ${base%$ext}tsv`
		
	if [ ! -f "$tsv" ]; then 
		echo could not find tsv file for $fn_bam 
		echo in folder $open_bis
		echo tsv: ${base%$ext}tsv
		exit 1
	fi
	
	echo $tsv

}


function merge_samples_tsv()
{

	if [ -z $1 -o -z $2 -o -z $3 -o -z $4 ]; then
		echo "Usage: $0 <bamDir> <file_extension> <fastqDir> <outputDir> <statusDir> [ --run-opt <submit> --no-braces <true/false> --sample-name <samples> --rmduplicates <true/false>]"
	fi
	
	indir=$1
	ext=$2
	open_bis=$3
	outdir=$4
	statusDir=$5
	
	for i in {1..5}
	do
		shift;
	done
	
	run_opt=local
	no_braces=false
	samples=""
	rmdup=false
	
	while [ 1 -lt $# ]; do
		if [ $1 == "--run-opt" ]; then
			shift;
			run_opt=$1;
		elif [ $1 == "--no-braces" ]; then
			shift;
			no_braces=$1;
		elif [ $1 == "--sample-name" ]; then
			shift;
			samples=$1;
		elif [ $1 == "--rmduplicates" ]; then
			shift;
			rmdup=$1;
		else
			echo did not understand arg $1
			shift
		fi;
	done
	
	
	mkdir -p $outdir
	
	if [ -z $samples ]; then
		for fn_bam in $indir/*$ext; do
		
			tsv=`get_tsv_for_bam $fn_bam $open_bis $ext`

			echo $fn_bam
			echo $tsv
			local sname=`extract $tsv EXTERNAL_SAMPLE_NAME`
			samples="$samples $sname"
			
		done
	fi
	
	if $no_braces; then
		samples=`echo $samples | sed 's/(//g' | sed 's/)//g' | sort -u`
	fi
	
	
	for sample_name in `echo $samples | sed 's/ /\n/g' | sort -u `; do 
		
		fn_result=$outdir/${sample_name}_mm.bam
	
		if [ -f $fn_result ]; then 
			echo "bam file $fn_result already exists."
			continue; 
		fi
	
		if [ $run_opt == "submit" ]; then
			echo "submit"
			log_dir=$outdir/logs/
			dir_=`pwd`
			mkdir -p $log_dir
			fn_out=${log_dir}merge_${sample_name}.o
			fn_err=${log_dir}merge_${sample_name}.e
			echo "" > $fn_err
			echo "" > $fn_out
			num_jobs=1000
			while [ "$num_jobs" -gt "5" ]
			do
				sleep 30
				num_jobs=`qstat | grep merge | wc -l`
				echo num_jobs: $num_jobs
			done
			queue="-q regular.q"
			echo "cd $dir_ ; echo ${sample_name} > $fn_out; hostname >> $fn_out; $0 $indir $ext $open_bis $outdir $statusDir --no-braces $no_braces --sample-name $sample_name --rmduplicates $rmdup " | qsub $queue -o $fn_out -e $fn_err -N merge_${sample_name} -cwd -V
		else
	
			bams=""
			joint_header=$outdir/joint_header_${sample_name}.txt
			if [ -f $joint_header ]; then
				rm $joint_header
			fi
			
			for fn_bam in $indir/*$ext
			do

				tsv=`get_tsv_for_bam $fn_bam $open_bis $ext`
				name=`extract $tsv EXTERNAL_SAMPLE_NAME`
				if $no_braces; then
				     name=`echo $name | sed 's/(//g' | sed 's/)//g'`
				fi
				if [ "$sample_name" == "$name" ]; then
	
					if [ ! -f $joint_header ]; then
						$samtools view -H $fn_bam > $joint_header
					else
						$samtools view -H $fn_bam | grep ^@RG >> $joint_header
					fi
					bams="$bams $fn_bam"
				fi
			done
			cat $joint_header | uniq > ${joint_header%.txt}.unique.txt
			mv ${joint_header%.txt}.unique.txt $joint_header
	
			echo $sample_name
			echo $bams
			fn_result=$outdir/${sample_name}_mm.bam
			fn_log=$outdir/$sample_name.log
		
			echo merge bam files > $fn_log
			for f in $bams; do
				echo $f | tee -a $fn_log
				samtools view $f | wc -l | tee -a $fn_log
			done
	
			echo "samtools merge -h $joint_header $fn_result $bams"
			$samtools merge -h $joint_header $fn_result $bams

			ret=$?
			if [ "$ret" != "0" ]; then 
				echo samtools merge failed with return value $ret
				exit 0 
			fi
			
			echo $fn_result | tee -a $fn_log
			samtools view $fn_result | wc -l | tee -a $fn_log 
	
			$samtools view -bh -F 256 $fn_result > ${fn_result%_mm.bam}.bam
	
			if $rmdup; then
				$picard_tools/MarkDuplicates.jar	INPUT=${fn_result%_mm.bam}.bam \
									OUTPUT=${fn_result%_mm.bam}.rmdup.bam \
									METRICS_FILE=${fn_result%_mm.bam}.rmdup.bam_metrics.txt \
									REMOVE_DUPLICATES=true \
									CREATE_INDEX=true \
									QUIET=true
			fi
		fi
	done
}
