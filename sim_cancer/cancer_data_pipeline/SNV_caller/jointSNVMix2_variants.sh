#/bin/bash -l

function jointSNVMix2_variants()
{	
	if [ -z "$1" -o -z "$3" -o -z "$4" -o -z "$5" ]; then # $2 is allowed to be empty
	        echo "usage jointSNVMix2_variants <prefix> <ext> <fn_genome_fasta> <output_dir> <normal_bam> <bed_file> [ --min_base_qual <0> --min_map_qual <0> ]"
	        exit -1;
	fi

	local prefix=$1; shift
	local ext=$1; shift
	local fn_genome=$1; shift
	local out_dir=$1; shift
	local bam_norm=$1; shift
	local bed_file=$1; shift
	
	mkdir -p $out_dir

	local run_with_options=false
	local jsm_options="" # for name of the output file
	local min_base_qual=0
	local min_map_qual=0

	while [ 0 -lt $# ]; do
		if [ $1 == "--min_base_qual" ]; then
			shift;
			min_base_qual=$1
			run_with_options=true
			jsm_options="${jsm_options}_min_base_qual_${min_base_qual}"
		elif [ $1 == "--min_map_qual" ]; then
			shift;
			min_map_qual=$1
			run_with_options=true
			jsm_options="${jsm_options}_min_map_qual_${min_map_qual}"
		fi
		shift
	done

	echo "Running jsm with the following parameters:"
	echo "prefix = $prefix"
	echo "ext = $ext"
	echo "fn_genome = $fn_genome"
	echo "out_dir = $out_dir"
	echo "bam_norm = $bam_norm"
	echo "bed_file = $bed_file"
	echo "run_with_options = $run_with_options"
	echo "jsm_options = $jsm_options"
	echo "min_base_qual = $min_base_qual"
	echo "min_map_qual = $min_map_qual"

	if [ -f $fn_genome -a ! -f $fn_genome.fai ]; then
	        echo did not find index file $fn_genome.fai
	        echo $samtools faidx $fn_genome
	        $samtools faidx $fn_genome
	fi
	
	if [ -f "$prefix" ]; then
	        local files=$prefix
	else
	        local files="$prefix*$ext"
	fi
	
	for fn_bam in $files
	do
		local outParams=$out_dir/`basename $fn_bam | sed 's/.bam//'`${jsm_options}.jointSNVMix2_Params.txt
		local outSNVsRaw=$out_dir/`basename $fn_bam | sed 's/.bam//'`${jsm_options}.jointSNVMix2_SNVs_Raw.txt
		
		## train parameters
		if [ ! -f $outParams ]; then
			opts=" --min_base_qual $min_base_qual --min_map_qual $min_map_qual "
			opts="$opts --skip_size 150 "
			opts="$opts --model beta_binomial "
			#opts="$opts --positions_file $bed_file "
			opts="$opts --priors_file $pathToJointSNVMix2/config/beta_binomial.priors.cfg " 
			opts="$opts --initial_parameters_file $pathToJointSNVMix2/config/beta_binomial.params.cfg "
			opts="$opts $fn_genome $bam_norm $fn_bam "
			opts="$opts $outParams "

			echo "$jointSNVMix train $opts"
			$jointSNVMix train $opts
		else
			echo "$outParams already exists"
		fi

		## perform variant calling
		if [ ! -f $outSNVsRaw ]; then
			opts=" --min_base_qual $min_base_qual --min_map_qual $min_map_qual "
			opts="$opts --model beta_binomial "
			opts="$opts --out_file $outSNVsRaw "
			opts="$opts --parameters_file $outParams "
			##opts="$opts --somatic_threshold 0.1 "
			opts="$opts --post_process "
			opts="$opts $fn_genome $bam_norm $fn_bam "

			echo "$jointSNVMix classify $opts"
			$jointSNVMix classify $opts
		else
			echo "$outSNVsRaw already exists"
		fi

		## convert to vcf format
		fn_vcf=${outSNVsRaw%.txt}.vcf
		if [ ! -f $fn_vcf ] ; then
			`find . -name JointSNVMixToVcf` $outSNVsRaw $fn_vcf
		else
			echo JointSNVMix vcf file exists: $fn_vcf
		fi

	done
}
