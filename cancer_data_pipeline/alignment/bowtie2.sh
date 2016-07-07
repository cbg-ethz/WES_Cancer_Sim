
function run_botie2()
{
	echo run_bowtie2

	# parse input
	################################################################################	
	local ref=$1; shift
	local prot=$1; shift
	local fastq1=$1; shift
	local fastq2=$1; shift
	local outdir=$1; shift
	local out_fname=$1; shift # e.g. tsv file name

	# optional arguments
	local sample_name="sname"
	local library="000"
	local k=20
	local sensitivity="--very-sensitive"
	local lane="xxx"
	local flowcellID="yyy"
	local opts=""

	while [ 0 -lt $# ]; do
		if [ $1 == "--sample" ]; then
			shift;
			sample_name=$1
		elif [ $1 == "--k" ]; then
			shift;
			k=$1
		elif [ $1 == "--sensitivity" ]; then
			shift;
			sensitivity=$1
		elif [ $1 == "--lane" ]; then
			shift;
			lane=$1
		elif [ $1 == "--flowcellID" ]; then
			shift;
			flowcellID=$1
		else
			opts="$opts $1"
		fi;
		shift;
	done



	# create directories
	################################################################################	
	mkdir -p $outdir/logs/
	mkdir -p $outdir/bam/
	mkdir -p $outdir/status/

	# define output files
	################################################################################	
	fn_status=$outdir/status/`basename ${out_fname%.tsv}`.done 
	fn_log=$outdir/logs/`basename ${out_fname%.tsv}`.log 
	fn_bam=$outdir/bam/`basename ${out_fname%.tsv}`.bam

	if [ -f $fn_status ]; then
		echo status file $fn_status exists
		return 0
	fi
	if [ -f $fn_bam ]; then
		echo "bam file $fn_bam already exists"
		return 0
	fi

	# define bowtie alignment parameters and count input reads
	################################################################################	
	if [ $prot == "PAIRED_END" ]; then 
		fastq="-1 $fastq1 -2 $fastq2" 
		cnt1=`zcat $fastq1 | wc -l`
		cnt2=`zcat $fastq2 | wc -l`
		cnt=`echo "($cnt1 + $cnt2)/4" | bc`
		echo number of reads in both fastq files `echo "$cnt/(1000^2)" | bc`M
	else
		fastq="-U $fastq1"
		cnt1=`zcat $fastq1 | wc -l`
		cnt=`echo "$cnt1/4" | bc`
		echo number of reads in fastq file `echo "$cnt/(1000^2)" | bc`M
	fi

	# compute index if needed
	################################################################################	
	index_base=${ref}_bowtie_idx/
	if [ ! -d $index_base ]; then 
		mkdir -p $index_base
		index_prefix=idx
		echo $bowtie_index_tool $ref $index_base/$index_prefix
		time $bowtie_index_tool $ref $index_base/$index_prefix
	else
		echo index exists: $index_base
		tmp=`basename $index_base*.1.bt2`
		index_prefix=${tmp%.1.bt2}
	fi
	index_base=$index_base/$index_prefix

	# compute number of free CPUs
	nr_threads=`compute_capacity`

	# set options
	opts="$opts -q" # input is fastq file
	opts="$opts --phred33"	#scaling of phred score
	opts="$opts $sensitivity" # alignment parameter set
	opts="$opts -p $nr_threads" # number of cores
	opts="$opts -k $k" # report up to k alignments per read
	#opts="$opts --trim3 5" # trim fixed amout of NTs 

	# read group settings
	rg="--rg-id ${flowcellID}.${lane}"
	rg="$rg --rg LB:$library"
	rg="$rg --rg PL:illumina"
	re="$rg --rg PU:${flowcellID}.${lane}"
	#rg="$rg --rg SM:`echo $sample_name  | cut -f 2 -d '.'`"
	rg="$rg --rg SM:$sample_name"

	# run tool
	################################################################################	

	echo "$bowtie $opts $rg -x $index_base $fastq 2> $fn_log | $samtools view -bhS \"-\" >> $fn_bam "
	time $bowtie $opts $rg -x $index_base $fastq 2> $fn_log | $samtools view -bhS "-" >> $fn_bam 


	ret=$?
	if [ ! "$ret" == "0" ]; then 
		echo bowtie finished with error $ret
		exit -1
	fi

	# check counts
	################################################################################	
	echo "$samtools view $fn_bam | wc -l"
	cnt_bam=`$samtools view $fn_bam | wc -l`
	echo "found $cnt_bam alignments from $cnt reads" 

	# assert that at least 80% are aligned 
	# in case of bowtie all unaligned reads are also represented in the file
	# therefore, the counts should be identical
	cnt_tol=`echo "$cnt*0.8/1" | bc`

	if [ "$cnt_bam" -lt "$cnt_tol" ]; then
		echo alignment failed:
		echo expected at least $cnt_to alignments from $cnt reads
		echo exit...
		exit -1
	fi


	# create status file 
	################################################################################	
	echo "$fn_bam" > $fn_status
	echo ";$cnt_bam; alignments from ;$cnt; reads" >> $fn_status &

	# sort bam file
	################################################################################	
										## if we sort with picard tools, we dont have problems later when using picard tools (i.e. in merge)
	echo "$picard_tools/SortSam.jar INPUT=$fn_bam OUTPUT=${fn_bam%.bam}.sorted SORT_ORDER=coordinate"
	$picard_tools/SortSam.jar INPUT=$fn_bam OUTPUT=${fn_bam%.bam}.sorted.bam SORT_ORDER=coordinate
	
	return 0
}

