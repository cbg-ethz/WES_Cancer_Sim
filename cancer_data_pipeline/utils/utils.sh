#!/bin/bash

function check_gzip()
{
	fn_stats=$1
	shift

	for f in "$@"
	do
		gzip --test $f
		if [ $? == "1" ]; then 
			echo $f | tee -a $fn_stats
			echo test for proper gzip format failed | tee -a $fn_stats
			mv $f $f.fail
		fi
	done
}

function count_fastq()
{
	fn_stats=$1
	shift

	for f in "$@"
	do

		cnt1=`zcat $f | wc -l`
		cnt=`echo "$cnt1/4" | bc`
		echo number of reads in $f | tee -a $fn_stats 
		echo "$cnt (`echo "$cnt/(1000^2)" | bc`M)"
		echo $cnt >> $fn_stats
	done
}

function extract()
{
	local tsv=$1
	local key=$2
	local value=`cat $tsv | grep $key | head -n 1 | cut -f 2-10 | sed 's/ /_/g' | sed 's/(//' | sed 's/)//'`
	echo $value
}

function check_bam()
{
	local fn_bam=$1
	x=`samtools view $fn_bam 2>&1 | head -n 2 | grep "EOF marker is absent."`; 
	if [ -z "$x" ]; then 
		echo ok
	else 
		echo truncated
	fi
}

function compute_capacity()
{
	local maximum=$1
	if [ -z "$maximum" ]; then
		maximum=100
	fi

	local load=`uptime | cut -f 4 -d , | cut -f 2 -d :`
	local nproc=`cat /proc/cpuinfo | grep processor| wc -l`
	if [ -z $load ]; then load=$nproc ; fi
	local NSLOTS=`echo "if ($nproc - $load<2) 1 else {nslots=1+($nproc - $load+0.9)/1; nslots}" | bc`

	if [ ! -z "$maximum" -a "$NSLOTS" -gt "$maximum" ]; then 
		NSLOTS=$maximum
	fi
	echo $NSLOTS
}

function compute_stats()
{
	local openbis_dir=$1
	# count number of files
	local rna_single_end=0
	local rna_paired_end=0
	local dna_single_end=0
	local dna_paired_end=0
	for tsv in `find $openbis_dir -name \*.tsv`; do
		local app=`extract $tsv SEQUENCING_APPLICATION`
		local protocol=`extract $tsv END_TYPE`
	
		if [ "$app" == "RNA_SEQ" ]; then
			if [ $protocol == "PAIRED_END" ]; then 
				rna_paired_end=$(($rna_paired_end + 1))
			else
				rna_single_end=$(($rna_single_end + 1))
			fi
		fi
	
		if [ "$app" == "GENOMIC_DNA_SEQ" ]; then
			if [ $protocol == "PAIRED_END" ]; then 
				dna_paired_end=$(($dna_paired_end + 1))
			else
				dna_single_end=$(($dna_single_end + 1))
			fi
		fi
	done
	echo
	echo DNA-Seq:
	echo $dna_single_end single end and $dna_paired_end paired end experiments
	echo
	echo RNA-Seq:
	echo $rna_single_end single end and $rna_paired_end paired end experiments
	echo
}

function run_fastqc()
{
	local fastq_dir=$1
	local out=$2
	local n_threads=$3

	mkdir -p $out
	for f in $fastq_dir/*.fastq.gz $fastq_dir/*.fastq ; do 

		
		if [ ! -f $f ]; then 
			echo "$f does not exist."
			continue
		fi

		echo run_fastqc: $fastqc on file $f
		local bn=`basename ${f%.fastq.gz}`
		local summary=$out/${bn%.fastq}_fastqc/summary.txt
		if [ ! -f $summary ]; then 
			if [ -z $n_threads ]; then
				echo $fastqc $f -o $out
				$fastqc $f -o $out
			else
				echo $fastqc $f -o $out -t $n_threads
                                $fastqc $f -o $out -t $n_threads
			fi
		fi
		cat $summary
	done
}

