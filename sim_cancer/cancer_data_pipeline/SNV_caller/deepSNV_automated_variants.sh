#!/bin/bash -l

function deepSNV_automated_variants()
{
	if [ -z "$1" -o -z "$2" -o -z "$3" -o -z "$4" ]; then #
		echo "usage deepSNV_automated_variants <bed_file> <tumor_bam> <normal_bam> <output_dir> [ --also-vcf <true/false> --run-opt <submit/locally> --minBaseQ <25> --estimateDispersion --overdispersion <100> --estimateDirichlet --alternative <two.sided> ]"
		exit -1;
	fi

	local bed_file=$1; shift
	local tumor_bam=$1; shift
	local normal_bam=$1; shift
	local out_dir=$1; shift
		
	local real_bed_file=${bed_file%.bed}_woFirstTwoLines_woLastColumn.bed
	
	local also_vcf=true
	local minBaseQ=25 # default for deepSNV
	local estimateDispersion=false
	local overdispersion=100 # default for deepSNV model betabinomial
	local estimateDirichlet=false
	local alternative="two.sided" 
	
	opts=""
	parameter_setting="alternative_${alternative}" # get one string that summarizes all the parameters set - and only specifically mentions them if they are non-default values, except for alternative

	while [ 0 -lt $# ]; do
		if [ $1 == "--also-vcf" ]; then
			shift;
			also_vcf=$1
		elif [ $1 == "--minBaseQ" ]; then
			shift;
			minBaseQ=$1
			if [ $minBaseQ != 25 ]; then
				parameter_setting="${parameter_setting}_minBaseQ${minBaseQ}"
				opts="$opts --minBaseQ $minBaseQ "
			fi
		elif [ $1 == "--estimateDispersion" ]; then
			estimateDispersion=true
			opts="$opts --estimateDispersion "
			parameter_setting="${parameter_setting}_estDispersionTRUE"
		elif [ $1 == "--overdispersion" ]; then
			shift;
			overdispersion=$1
			if [ $overdispersion != 100 ]; then
				parameter_setting="${parameter_setting}_overdisp${overdispersion}"
				opts="$opts --overdispersion $overdispersion "
			fi
		elif [ $1 == "--estimateDirichlet" ]; then
			estimateDirichlet=true
			opts="$opts --estimateDirichlet "
			parameter_setting="${parameter_setting}_estDirichletTRUE"
		elif [ $1 == "--alternative" ]; then
			shift;
			alternative=$1;
			opts="$opts --alternative $alternative "
		fi;
		shift;
	done

	echo "run deepSNV..."
	echo "opts = $opts "
	
	mkdir -p $out_dir
	
	if [ $((`cat $bed_file | head -n 100 | wc -l` *3)) -eq `cat $bed_file | head -n 100 | wc -w` ]; then 
		# this bed file has exactly three columns 
		real_bed_file=$bed_file
	else
		# assume this is a bed file from agilent with two header lines and 4 columns
		if [ ! -f $real_bed_file  ]; then
			#cat $bed_file | awk 'NR > 2' | awk '{$NF=""; print $0}'  > $real_bed_file ## I leave this out because if this is invoked multiple times, deepSNV might start on a wrong bed file
			echo "Error: The bed file has header lines, which deepSNV will not accept. Please generate the bed file with no header lines and just the regions."
			exit 0
		fi
	fi
	 
	if [ ! -f $tumor_bam.bai ]; then 
		echo "$samtools index $tumor_bam"
		$samtools index $tumor_bam
	fi
	if [ ! -f $normal_bam.bai ]; then 
		echo "$samtools index $normal_bam"
		$samtools index $normal_bam
	fi
	
	deepSNV_script=`find $SNV_DIR -name deepSNV_automated.R`
	 
	command="$deepSNV_script $real_bed_file $tumor_bam $normal_bam $out_dir $opts"
	echo $command
	eval $command
	
	if $also_vcf; then
		output=$out_dir/`basename ${tumor_bam%.bam}`_`basename ${normal_bam%.bam}`_${parameter_setting}.deepSNV.txt
		fn_vcf=${output%.txt}.vcf
		if [ ! -f $fn_vcf ] ; then 
			`find $SNV_DIR -name deepSNVToVcf` $output $fn_vcf 
		else
			echo deepSNV vcf file exists: $fn_vcf
		fi
	fi
}
