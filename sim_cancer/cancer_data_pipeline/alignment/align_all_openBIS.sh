#!/bin/bash

if [ -z $1 ]; then 
	echo "Usage: $0 <input_dir> <output_dir> [local]"
fi

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`
source `find ../../ -name paths.sh`
source `find ${dir_} -name paths.sh`
source `find ${dir_} -name utils.sh`
source `find ${dir_} -name genome.sh`
source `find ${dir_} -name bowtie2.sh`
source `find ${dir_} -name submit.sh`

openbis_dir=$1
outdir=$2
run_opt=$3

mkdir -p $outdir

################################################################################
# if arg3 not specified, then submit all jobs to cluster
################################################################################
if [ -z "$run_opt" -o "$run_opt" == "submit" ]; then 
	echo submit $run_opt 
	cnt=0
	log_dir=$outdir/logs/
	dir_=`pwd`
	mkdir -p $log_dir
	for f in `find $openbis_dir -name \*.tsv`; do
		cnt=$((cnt + 1))
		fn_out=${log_dir}align_${cnt}_`basename $f`.o
		command="cd $dir_ ; echo $f > $fn_out; hostname >> $fn_out; $0 `dirname $f` $outdir ToGrid"
		submit "$command" --fn-out $fn_out --tag align_tsv$cnt
		wait_for_jobs align_tsv --max-num 15
	done
	wait_for_jobs align_tsv
	exit 0
fi
################################################################################


################################################################################
# count number of jobs 
################################################################################
compute_stats $openbis_dir
################################################################################

if [ $run_opt == "single" ]; then 
	tsv_files=$openbis_dir/*.tsv
else
	tsv_files=`find $openbis_dir -name \*.tsv`
fi

for tsv in $tsv_files; do
	echo $tsv

	app=`extract $tsv SEQUENCING_APPLICATION`
	protocol=`extract $tsv END_TYPE`
	sample_name=`extract $tsv EXTERNAL_SAMPLE_NAME`
	readlen=`extract $tsv CYCLES_REQUESTED_BY_CUSTOMER`
	organims=`extract $tsv NCBI_ORGANISM_TAXONOMY`
	#lane=`extract $tsv LANECOUNT`
	lane=$(echo `basename $tsv` | tr '_' ' ' | awk '{print $10}' | sed 's/L00//')
	library=$(echo `basename $tsv` | tr '_' ' ' | awk '{print $1"_"$2"_"$3}')
	flowcellID=$(echo `basename $tsv` | tr '_' ' ' | awk '{print $7}')

	# get genome using function in genome.sh
	ref=`get_genome $organims`
	if [ "$ref" == "unknown" ]; then 
		echo did not recognice organism: $organism, continue... 
		continue
	else
		echo using genome $ref
	fi
	
	dd=`dirname $tsv`
	echo $dd
	for f in `ls $dd | grep _R1`; do 
		fastq1=$dd/$f
		if [ -z $fastq1 ]; then
			num_fastq=`ls $dd/*fastq.gz | wc -l`
			if [ $num_fastq == "1" ]; then 
				fastq1=$dd/*fastq.gz
			fi
		fi
		if [ -z $fastq1 ]; then
			echo $app
			echo did not find fastq file in $dd
			exit -1
		fi
		fastq2=`echo $fastq1 | sed 's/_R1/_R2/'`
		opts=""

		if [ ! "$app" == "GENOMIC_DNA_SEQ" ]; then
			if [ $protocol == "PAIRED_END" ]; then 
				echo RNA_Seq paired, continue
				continue
				#echo run_palmapper $ref $protocol $fastq1 $fastq2 $sample_name $outdir $lane $tsv
			else
				echo RNA_Seq single, continue
				continue
				#echo single end: $tsv
				#continue
				#echo run_palmapper $ref $protocol $fastq1 "xxx" $sample_name $outdir $lane $tsv
			fi
		else
			if [ $protocol == "PAIRED_END" ]; then 
				if [ "$fastq2" == "$fastq1" ]; then 
					echo "could not find file: $fastq2"
					exit 0
				fi
				run_botie2 $ref $protocol $fastq1 $fastq2 $sample_name $outdir $library $tsv "$opts" $lane $flowcellID
			else
				run_botie2 $ref $protocol $fastq1 "xxx" $sample_name $outdir $library $tsv "$opts" $lane $flowcellID
			fi
		fi
	done
	if [[ -z `ls $dd | grep _R1` &&  ! -z `ls $dd | grep _R2` ]]; then ## when it was clipped before, we have the unpaired R1 and R2 reads in separate folders; without this additional if condition, the unpaired R2 file would not be aligned
		for f in `ls $dd | grep _R2`; do
			fastq1=$dd/$f
	                if [ -z $fastq1 ]; then
	                        echo $app
	                        echo did not find fastq file in $dd
	                        exit -1
	                fi
			if [ "$app" == "GENOMIC_DNA_SEQ" ]; then
				if [ ! $protocol == "PAIRED_END" ]; then
					run_botie2 $ref $protocol $fastq1 "xxx" $sample_name $outdir $library $tsv "$opts" $lane $flowcellID
				fi
			fi
 
		done
	fi

done

