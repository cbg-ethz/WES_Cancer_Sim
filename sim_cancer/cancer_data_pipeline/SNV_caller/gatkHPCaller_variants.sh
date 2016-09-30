#/bin/bash -l

function gatkHPCaller_variants()
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
	
	if [ -z "$ext" -a -f "$prefix" ]; then
		local files=$prefix
	else
		local files="$prefix*$ext"
	fi
	
	for fn_bam in $files ; do
		local outputTU=$out_dir/`basename $fn_bam | sed 's/.bam//'`.gatkHPCaller_SNVs.raw.vcf
		local outputNO=$out_dir/`basename $fn_bam | sed 's/.bam//'`.NOsample.gatkHPCaller_SNVs.raw.vcf	
		local outputpaired=$out_dir/`basename $fn_bam | sed 's/.bam//'`.paired_TU_NO.gatkHPCaller_SNVs.raw.vcf

		if [ ! -f $outputTU ]; then

			# by default read with mapping qualitiy of 255 are filtered out by gatk. However, bowtie assigns mapping quality of 
			# 255 to all reads with only a single alignment. With this option we modify this mapping quality
			rq_opt="-rf ReassignOneMappingQuality -RMQF 255 -RMQT 254"
			#$GATK		-R $fn_genome -T HaplotypeCaller -I $fn_bam -I $bam_norm -o $output $rq_opt
			command="$GATK -R $fn_genome -T HaplotypeCaller -I $fn_bam -o $outputTU $rq_opt"
			echo $command
			eval $command


			command="$GATK -R $fn_genome -T HaplotypeCaller -I $bam_norm -o $outputNO $rq_opt"
			echo $command
			eval $command
		else
			echo "gatk output file already exists $outputTU."
		fi 	

		if [ ! -f $outputpaired ]; then

			# by default read with mapping qualitiy of 255 are filtered out by gatk. However, bowtie assigns mapping quality of 
			# 255 to all reads with only a single alignment. With this option we modify this mapping quality
			rq_opt="-rf ReassignOneMappingQuality -RMQF 255 -RMQT 254"
			command="$GATK -R $fn_genome -T HaplotypeCaller -I $fn_bam -I $bam_norm  -o $outputpaired $rq_opt"
			echo $command
			eval $command
		else
			echo "gatk output file already exists $outputpaired."
		fi 	

	done
}
