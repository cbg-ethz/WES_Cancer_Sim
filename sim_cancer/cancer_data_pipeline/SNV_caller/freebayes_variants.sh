#!/bin/bash

function freebayes_variants()
{
	local prefix=$1; shift; 
	local ext=$1; shift; 
	local fn_genome=$1; shift; 
	local out_dir=$1; shift;
	
	local min_alternate_count=""
	local min_alternate_fraction=""
	local pooled_continuous=""
	local min_mapping_quality=""
	local threads=""
	local freebayes_options=""

	while [ 0 -lt $# ]; do
		if [ $1 == "--min-alternate-count" ]; then
			shift;
			min_alternate_count=$1
			freebayes_options="${freebayes_options}_min_alternate_count_${min_alternate_count}"
		elif [ $1 == "--min-alternate-fraction" ]; then
			shift;
			min_alternate_fraction=$1;
			freebayes_options="${freebayes_options}_min_alternate_fraction_${min_alternate_fraction}"
		elif [ $1 == "--pooled-continuous" ]; then
			shift;
			pooled_continuous=$1
			freebayes_options="${freebayes_options}_pooled_continuous"
		elif [ $1 == "--min-mapping-quality" ]; then
			shift;
			min_mapping_quality=$1
			freebayes_options="${freebayes_options}_min_mapping_quality_${min_mapping_quality}"
		elif [ $1 == "--freebayes-parallel" ]; then
			shift;
			threads=$1
		else
			shift;
		fi;
	done
	
	if [ -z "$prefix" -o -z "$fn_genome" -o -z "$out_dir" ]; then 
		echo "usage freebayes_variants <prefix> <ext> <fn_genome_fasta> <output_dir> [ --min-alternate-count <min_alternate_count> --min-alternate-fraction <min_alternate_fraction> --pooled-continuous <pooled_continuous> --min-mapping-quality <min_mapping_quality> --freebayes-parallel <8> ]"
		exit -1;
	fi
	
	mkdir -p $out_dir
	
	if [ -f "$prefix" ]; then
		files=$prefix
	else
		files="$prefix*$ext"
	fi

	# to check arguments
	echo "These were the arguments given to freebayes_variants:"
	echo "prefix = $prefix"
	echo "ext = $ext"
	echo "fn_genome = $fn_genome"
	echo "out_dir = $out_dir"
	echo "min_alternate_count = $min_alternate_count"
	echo "min_alternate_fraction = $min_alternate_fraction"
	echo "pooled_continuous = $pooled_continuous"
	echo "min_mapping_quality = $min_mapping_quality"
	echo "threads = $threads"
	echo "freebayes_options = $freebayes_options"
	
	for fn_bam in $files
	do
		fn_out=${out_dir}/`basename ${fn_bam%.bam}`${freebayes_options}.freebayes.vcf
		if [ ! -f $fn_out ]; then 
			opts=" "
			if [ ! -z "$min_alternate_count" ]; then
				opts="$opts --min-alternate-count $min_alternate_count"
			fi
			if [ ! -z "$min_alternate_fraction" ]; then
				opts="$opts --min-alternate-fraction $min_alternate_fraction"
			fi
			if [ ! -z "$pooled_continuous" ]; then
				opts="$opts --pooled-continuous"
			fi
			if [ ! -z "$min_mapping_quality" ]; then
				opts="$opts --min-mapping-quality $min_mapping_quality"
			fi

			if [ -z "$threads" ]; then
				echo "$freebayes -f $fn_genome $opts $fn_bam > $fn_out"
				time $freebayes  -f $fn_genome $opts $fn_bam > $fn_out
				ret=$?
				if [ ! "$ret" -eq "0" ]
				then 
					echo "freebayes returned error code $ret"
					exit $ret
				fi
			else
				echo "$freebayes_parallel <($fasta_generate_regions $fn_genome.fai 100000) $threads -f $fn_genome $opts $fn_bam > $fn_out"
				$freebayes_parallel <($fasta_generate_regions $fn_genome.fai 100000) $threads -f $fn_genome $opts $fn_bam > $fn_out
				ret=$?
                                if [ ! "$ret" -eq "0" ]
                                then
                                        echo "freebayes returned error code $ret"
                                        exit $ret
                                fi
			fi
		else
			echo $fn_out exists
		fi
	done
}
