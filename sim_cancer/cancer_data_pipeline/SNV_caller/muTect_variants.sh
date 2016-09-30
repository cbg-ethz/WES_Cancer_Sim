#/bin/bash -l

function muTect_variants()
{

	source `find ${gitDir} -name paths.sh`
	
	local prefix=$1
	local ext=$2
	local fn_genome=$3
	local out_dir=$4
	local bam_norm=$5
	
	if [ -z "$1" -o -z "$3" -o -z "$4" -o -z "$5" ]; then # $2 is allowed to be empty
		echo "usage $0 <prefix> <ext> <fn_genome_fasta> <output_dir> <normal_bam>"
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
		output=$out_dir/`basename $fn_bam | sed 's/.bam//'`.muTect_SNVs_Raw.txt
		outputCoverage=$out_dir/`basename $fn_bam | sed 's/.bam//'`.muTect_coverage.wig.txt
		outputVCF=$out_dir/`basename $fn_bam | sed 's/.bam//'`.muTect_SNVs_Raw.vcf
	
		echo "MuTect output: $output"
		echo "Check space on node:"
		df /tmp/
	
		if [ ! -s $output ]; then
		
			echo run mutect
			local opts=""
			$mutect --input_file:normal $bam_norm --input_file:tumor $fn_bam \
					-R $fn_genome --out $output -vcf $outputVCF \
					--coverage_file $outputCoverage $opts
	
		else
			echo "MuTect output file already exists."
		fi 
		
	done
}
