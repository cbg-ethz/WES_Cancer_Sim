#/bin/bash -l

function somaticSniper_variants()
{
	local prefix=$1
	local ext=$2
	local fn_genome=$3
	local out_dir=$4
	local bam_norm=$5
	
	if [ -z "$1" -o -z "$3" -o -z "$4" -o -z "$5" ]; then # $2 is allowed to be empty
	        echo "usage somaticSniper_variants <prefix> <ext> <fn_genome_fasta> <output_dir> <normal_bam>"
	        exit -1;
	fi
	
	mkdir -p $out_dir
	
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
		local vcf_output=$out_dir/`basename $fn_bam | sed 's/.bam//'`.somaticSniper_SNVs_Raw.vcf
	
		if [ ! -f $vcf_output ]; then
			$somaticSniper -F vcf        -f	$fn_genome \
	                                        	$fn_bam \
	                                        	$bam_norm \
	                                        	$vcf_output
	
		else
			echo "Somatic sniper vcf-output file already exists."
		fi
	
	done
}
