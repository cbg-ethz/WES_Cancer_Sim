
function sinvict_variants()
{
	source `find ${gitDir} -name paths.sh`

	if [ -z "$1" -o -z "$2" -o -z "$3" -o -z "$4" -o -z "$5" ]; then
		echo "usage $0 <bamFile> <fn_genome_fasta> <output_dir> <normal_bam> <bed_file>"
		exit -1;
	fi
	
	bamFile=$1; shift;
	fn_genome=$1; shift;
	out_dir=$1; shift;
	bam_norm=$1;shift;
	bed_file=$1; shift

	sinvict_options=""
	sinvict_opts_short=""

	while [ 0 -lt $# ]; do
		sinvict_options=" $sinvict_options $1 "
		shift
	done

	if [ ! -z "$sinvict_options" ]; then
		echo "Running sinvict with the options: $sinvict_options"
		sinvict_opts_short=`echo $sinvict_options | sed 's/-//g' | sed 's/ //g'`
	fi

	mkdir -p $out_dir

	# first: bam-readcount to get the input files:
	readcountFolderTU=$out_dir/readCounts_`basename $bamFile`/
	mkdir -p $readcountFolderTU
	readcountout=$readcountFolderTU/`basename $bamFile`_readcount.txt
	if [ ! -f $readcountout ]; then
		command="$bam_readcount -w 1 -f $fn_genome $bamFile -l $bed_file > $readcountout"
		echo $command
		eval $command
	fi

	readcountFolderNO=$out_dir/readCounts_`basename $bam_norm`_${RANDOM}/
	mkdir -p $readcountFolderNO
	readcountout=$readcountFolderNO/`basename $bam_norm`_readcount.txt
	if [ ! -f $readcountout ]; then
	        command="$bam_readcount -w 1 -f $fn_genome $bam_norm -l $bed_file > $readcountout"
	        echo $command
	        eval $command
	fi

	sinvictOutTU=$out_dir/sinvict${sinvict_opts_short}_`basename $bamFile`/
	if [ ! -f $sinvictOutTU/calls_level1.sinvict ]; then
		mkdir -p $sinvictOutTU
		command="$sinvict $sinvict_options  -t $readcountFolderTU/ -o $sinvictOutTU"
		echo $command
		eval $command
	fi

	sinvictOutNO=$out_dir/sinvict${sinvict_opts_short}_`basename $bam_norm`_${RANDOM}/
	if [ ! -f $sinvictOutNO/calls_level1.sinvict ]; then
		mkdir -p $sinvictOutNO
		command="$sinvict $sinvict_options  -t $readcountFolderNO/ -o $sinvictOutNO"
		echo $command
		eval $command
	fi

	# subtract all variants that were found in the normal, and rewrite into vcf format:
	sinvict2somaticVCF=`find ${gitDir} -name sinvict2somaticVCF.py`
	tumorSiNCIVTfilesOrdered=""
	for i in `seq 1 6`; do
		tumorSiNCIVTfilesOrdered="$tumorSiNCIVTfilesOrdered  $sinvictOutTU/calls_level${i}.sinvict "
	done
	outSomaticVCF=$out_dir/sinvict${sinvict_opts_short}_`basename ${bamFile%.bam}`_vs_`basename ${bam_norm%.bam}`_somatic.vcf
	command="$python $sinvict2somaticVCF $tumorSiNCIVTfilesOrdered $sinvictOutNO/calls_level1.sinvict $outSomaticVCF $fn_genome"
	if [ ! -f $outSomaticVCF ]; then
		echo $command
		eval $command
	fi	
}
