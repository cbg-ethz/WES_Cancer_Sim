#!/bin/bash

function samtools_1_2_variants()
{
	source `find ${gitDir} -name paths.sh`
	
	local prefix=$1; shift; 
	local ext=$1; shift; 
	local fn_genome=$1; shift; 
	local out_dir=$1; shift;
	local bam_norm=$1; shift;

	local sam_options=""
	local E_opt=" -E "
	local paired=false

	while [ 0 -lt $# ]; do
		if [ $1 == "--opts" ]; then
			shift
			sam_options=$1
		elif [ $1 == "--paired" ]; then
			paired=true
		fi
                shift
	done
	
	echo "sam_options = $sam_options"
	if [[ "$sam_options" == *"--no-E"* ]]; then
		E_opt=""
		sam_options=$(echo $sam_options | sed 's/--no-E/ /' )
	fi
	echo "paired = $paired"

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
		if [ ! -z "$sam_options" ]; then
			fn_outTU=${out_dir}/`basename ${fn_bam%bam}`option_`echo $sam_options | sed 's/-//g' | sed 's/ //g'`_samvar_1_2.vcf.gz
			fn_outNO=${out_dir}/`basename ${fn_bam%bam}`option_`echo $sam_options | sed 's/-//g' | sed 's/ //g'`_samvar_1_2.vcf.NOsample.gz
			fn_out_paired=${out_dir}/`basename ${fn_bam%bam}`option_`echo $sam_options | sed 's/-//g' | sed 's/ //g'`_samvar_1_2.vcf.paired_TU_NO.gz
		else
			fn_outTU=${out_dir}/`basename ${fn_bam%bam}`_samvar_1_2.vcf.gz
			fn_outNO=${out_dir}/`basename ${fn_bam%bam}`_samvar_1_2.vcf.NOsample.gz
			fn_out_paired=${out_dir}/`basename ${fn_bam%bam}`_samvar_1_2.vcf.paired_TU_NO.gz
		fi

		if [ -z "$E_opt" ]; then
			fn_outTU=${fn_outTU%_samvar_1_2.vcf.gz}_noE_samvar_1_2.vcf.gz
			fn_outNO=${fn_outTU%_samvar_1_2.vcf.gz}_noE_samvar_1_2.vcf.NOsample.gz
			fn_out_paired=${fn_outTU%_samvar_1_2.vcf.gz}_noE_samvar_1_2.vcf.paired_TU_NO.gz
		else
			fn_outTU=${fn_outTU%_samvar_1_2.vcf.gz}_withE_samvar_1_2.vcf.gz
			fn_outNO=${fn_outTU%_samvar_1_2.vcf.gz}_withE_samvar_1_2.vcf.NOsample.gz
			fn_out_paired=${fn_outTU%_samvar_1_2.vcf.gz}_withE_samvar_1_2.vcf.paired_TU_NO.gz
		fi

		if [ ! -f $fn_outTU ]; then 
			touch $fn_outTU
			# option -E results in more sensitive predictions
			opts=""
			##opts=" $opts -C50 " ##The default is 0
			if [[ ! $sam_options == *"-B"*  ]]; then #  The -B option cannot be combined with -E
				opts=" $opts $E_opt " ##Recalculate BAQ on the fly, ignore existing BQ tags (BAQ = Base Alignment Quality)
			fi
			opts=" $opts --output-tags DP4,DP "
			opts=" $opts $sam_options "
			command="$samtools_1_2 mpileup $opts -uf $fn_genome $fn_bam | $bcftools_1_2 call -mv -Oz > $fn_outTU"
			echo $command
			eval $command
			ret=$?
			if [ ! "$ret" -eq "0" ]
			then 
				echo "samtools_1_2 returned error code $ret"
				exit $ret
			fi
		else
			echo $fn_outTU exists
		fi

		if [ ! -f $fn_outNO ]; then
			command="$samtools_1_2 mpileup $opts -uf $fn_genome $bam_norm | $bcftools_1_2 call -mv -Oz > $fn_outNO"
			echo $command
			eval $command
			ret=$?
			if [ ! "$ret" -eq "0" ]
			then
				echo "samtools_1_2 returned error code $ret"
				exit $ret
			fi
		else
			echo $fn_outNO exists
		fi


		if $paired; then
			if [ ! -f $fn_out_paired ]; then
				# option -E results in more sensitive predictions
				opts=""
				##opts=" $opts -C50 " ##The default is 0
				if [[ ! $sam_options == *"-B"*  ]]; then #  The -B option cannot be combined with -E
					opts=" $opts $E_opt " ##Recalculate BAQ on the fly, ignore existing BQ tags (BAQ = Base Alignment Quality)
				fi
				opts=" $opts --output-tags DP4,DP "
				opts=" $opts $sam_options "
				command="$samtools_1_2 mpileup $opts -uf $fn_genome $bam_norm $fn_bam | $bcftools_1_2 call -mv -Oz > $fn_out_paired"
				echo $command
				eval $command
				ret=$?
				if [ ! "$ret" -eq "0" ]
				then 
					echo "samtools_1_2 returned error code $ret"
					exit $ret
				fi			
			else
				echo $fn_out_paired exists
			fi
		fi
	done
}
