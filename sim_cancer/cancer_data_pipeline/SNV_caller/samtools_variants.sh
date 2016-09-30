#!/bin/bash

function samtools_variants()
{
	local prefix=$1; shift; 
	local ext=$1; shift; 
	local fn_genome=$1; shift; 
	local out_dir=$1; shift;
	
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
		fn_out=${out_dir}/`basename ${fn_bam%bam}samvar.bcf`
		if [ ! -f $fn_out ]; then 
			# option -E results in more sensitive predictions
			opts=" "
			##opts=" $opts -C50 " ##The default is 0
			opts=" $opts -E " ##Recalculate BAQ on the fly, ignore existing BQ tags (BAQ = Base Alignment Quality)
			echo "$samtools mpileup $opts -uf $fn_genome $fn_bam | $bcftools view -bvcg - > $fn_out"
			time $samtools mpileup $opts -uf $fn_genome $fn_bam | $bcftools view -bvcg - > $fn_out
			ret=$?
			if [ ! "$ret" -eq "0" ]
			then 
				echo "samtools returned error code $ret"
				exit $ret
			fi
		else
			echo $fn_out exists
		fi
	done
	
}
