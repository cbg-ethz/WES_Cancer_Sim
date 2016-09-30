#/bin/bash -l

function varscan2_variants()
{
	source `find ${gitDir} -name paths.sh`

	local prefix=$1; shift
	local ext=$1; shift;
	local fn_genome=$1; shift;
	local out_dir=$1; shift;
	local bam_norm=$1; shift;

	if [ -z "$prefix" -o -z "$fn_genome" -o -z "$out_dir" -o -z "$bam_norm" ]; then # $2 is allowed to be empty
		echo "prefix: $prefix"
		echo "ext: $ext"
		echo "fn_genome: $fn_genome"
		echo "out_dir: $out_dir"
		echo "fn_bam_NO: $bam_norm"
		echo "usage varscan2_variants <prefix> <ext> <fn_genome_fasta> <output_dir> <normal_bam> [<varscan_opts>]"
		exit -1;
	fi

	local run_with_options=false
	local varscan_options="" # for name of the output file

	while [ 0 -lt $# ]; do
		run_with_options=true
		varscan_options="$varscan_options $1"
		shift
	done

	echo "Running varscan2 with the options: $varscan_options"
	
	mkdir -p $out_dir
	
	if [ -f $fn_genome -a ! -f $fn_genome.fai ]; then
		echo did not find index file $fn_genome.fai
		echo $samtools faidx $fn_genome
		$samtools faidx $fn_genome
	fi
	
	if [ -f "$prefix" ]; then
		files=$prefix
	else
		files="$prefix*$ext"
	fi
	
	for fn_bam in $files
	do
		bam_tumor_prefix=$out_dir/`basename $fn_bam`_`echo $varscan_options | sed 's/-//g' | sed 's/ //g'`
		bam_normal_prefix=$out_dir/`basename $bam_norm`_`echo $varscan_options | sed 's/-//g' | sed 's/ //g'`
		tumor_pileup=${bam_tumor_prefix%.bam}.$RANDOM.pileup
		normal_pileup=${bam_normal_prefix%.bam}.$RANDOM.pileup # make sure different processes don't iterfere 
	
	 	if [ ! -f ${bam_tumor_prefix%.bam}_varscan2.txt.snp.Somatic.hc ]; then

			if [ ! -f $normal_pileup ]; then
				echo "$samtools mpileup -f $fn_genome $bam_norm > $normal_pileup"
				$samtools mpileup -f $fn_genome $bam_norm > $normal_pileup
			fi
	 		if [ ! -f $tumor_pileup  ]; then
				opts=" "
				##opts=" $opts -Q 5 " ## default is 13
	 			echo "$samtools mpileup -f $opts $fn_genome $fn_bam > $tumor_pileup"
	 			$samtools mpileup -f $opts $fn_genome $fn_bam > $tumor_pileup
	 		fi
	 
			if [[ $varscan_options ==  *--min-tumor-freq* ]]; then # this is an option for the processSomatic command
				command="$varscan2 somatic $normal_pileup $tumor_pileup ${bam_tumor_prefix%.bam}_varscan2.txt"
				echo $command
				eval $command
				command="$varscan2 processSomatic ${bam_tumor_prefix%.bam}_varscan2.txt.snp $varscan_options"
				echo $command
				eval $command
			else
				command="$varscan2 somatic $normal_pileup $tumor_pileup ${bam_tumor_prefix%.bam}_varscan2.txt $varscan_options"
		 		echo $command 
		 		eval $command
				command="$varscan2 processSomatic ${bam_tumor_prefix%.bam}_varscan2.txt.snp"
				echo $command
				eval $command 
			fi
			## -strand-filter 1: removes variants with greater than 90% strand bias _______________________________________________________
	  		#__ processSomatic: separates the somatic output file by: Germline, Somatic, LOH. Somatic mutations will further be classified as high-confidence (.hc) or low-confidence (.lc)
	 		rm $tumor_pileup
	 		rm $normal_pileup
	 	else
	 		echo "${bam_tumor_prefix%.bam}_varscan2.txt.snp.Somatic.hc exists"
	 	fi
	
		if [ ! -f ${bam_tumor_prefix%.bam}_varscan2.txt.snp.Somatic.hc.vcf ]; then
			echo "`find ${gitDir} -name varscan2ToVcf.sh` ${bam_tumor_prefix%.bam}_varscan2.txt.snp.Somatic.hc $fn_genome"
			time `find ${gitDir} -name varscan2ToVcf.sh` ${bam_tumor_prefix%.bam}_varscan2.txt.snp.Somatic.hc $fn_genome
		fi
	
		if [ ! -f ${bam_tumor_prefix%.bam}_varscan2.txt.snp.Somatic.vcf ]; then
			echo "`find ${gitDir} -name varscan2ToVcf.sh` ${bam_tumor_prefix%.bam}_varscan2.txt.snp.Somatic $fn_genome"
			time `find ${gitDir} -name varscan2ToVcf.sh` ${bam_tumor_prefix%.bam}_varscan2.txt.snp.Somatic $fn_genome
		fi

	done
}
