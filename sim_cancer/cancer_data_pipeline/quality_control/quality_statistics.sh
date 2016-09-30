#!/bin/bash

# Description:
# This script performs quality checks after the pipeline is finished. For the bam files that are given as input, 
# it computes an output file which contains the following information: # reads in the original fastq(s); # reads
# that survived the clipping; % of reads that survived the clipping; # reads aligned; % of reads aligned; minimum,
# maximum, mean, and median coverage of the genes; % of bases in the genes that are not covered at all, covered with
# at least 20 reads; % of reads that overlap the genes.
# It also executes fastqc.

# Inputs:
# BamFiles = One or several bam files for which the quality control should be conducted. They are expected to be named <sample_name><ext>
# ext = the extension that comes in the name of the bam file after the sample name. For example ".bam" or ".rmdup.bam"
# stats_dir = During the clipping, a folder called "stats" was created which contains read counts for each pair of fastq files
# gtf_file_genes = the gtf file that contains only the gene boundaries (exons, introns, etc. removed) OR the BED file with the targeted regions of the exome kit
# output_dir = the output folder for all the quality control files
# chromInfo = the file (often called chromInfo.txt), which lists the names of all chromosomes together with the size of the chromosome in basepairs; 

# Outputs:
# Within the output folder there will be the following directories:
# coverage = for each sample the files used to determine the coverage; most of them generated with BEDtools
# <sample_name>_fastqc = the folder which serves as output from the program fastqc
# In the output folder there will be the following files:
# <sample_name>_read_statistics.txt = containing a summary of all the quality information


function quality_statistics()
{

        if [ -z "$1" -o -z "$2" -o -z "$3" -o -z "$4" -o -z "$5" -o -z "$6" ]; then
                echo "usage quality_statistics <BamFiles> <ext> <stats_dir_from_clip_trim> <gtf_file_genes> <output_dir> <chromInfo> <tsv_dir> [ --run-opt <submit/locally> ]"
                exit -1;
        fi

        local bamFileS=$1; shift;
	local ext=$1; shift;
	local stats_dir=$1; shift;
        local gtf_file=$1; shift;
        local out_dir=$1; shift;
	local chromInfo=$1; shift;
	local tsv_dir=$1; shift;

        local run_opt=locally

        local i=1
        while [ $i -lt $# ]; do
                if [ $1 == "--run-opt" ]; then
                        shift;
                        run_opt=$1;
                else
                        echo did not understand arg $1
                        shift;
                fi;
        done
        echo "Running quality_statistics using the following inputs:"
        echo "bamFileS = $bamFileS"
	echo "ext = $ext"
	echo "stats_dir = $stats_dir"
	echo "gtf_file_genes = $gtf_file"
        echo "out_dir = $out_dir"
	echo "chromInfo = $chromInfo"
	echo "tsv_dir = $tsv_dir"
        echo "Submit to grid or run locally = $run_opt"

        if [ "$run_opt" == "submit" ]; then
                echo submit $run_opt

                #####################
                # run it
                local log_dir=$out_dir/logs/
                local dir_=`pwd`
                mkdir -p $log_dir
                local cnt=0

                for fn_bam in `echo $bamFileS`; do
                        echo current bam $fn_bam
			sample_name=`basename ${fn_bam%$ext}`
			echo "sample_name = $sample_name"
	
			if [ -f $out_dir/${sample_name}_read_statistics.txt ]; then
				echo "$out_dir/${sample_name}_read_statistics.txt already exists"
				continue
			fi

                        cnt=$((cnt + 1))
                        local fn_out=${log_dir}quality_statistics${cnt}_`basename ${fn_bam%.bam}`.o
                        local fn_err=${log_dir}quality_statistics${cnt}_`basename ${fn_bam%.bam}`.e
                        echo "" > $fn_err
                        echo "" > $fn_out
                        num_jobs=1000
                        while [ "$num_jobs" -gt "5" ]
                        do
                                sleep 10
                                num_jobs=`qstat | grep quality | wc -l`
                                echo num_jobs: $num_jobs
                        done

			queue="-q regular.q"
                        #queue="-q regular.q@bs-dsvr02,regular.q@bs-dsvr04,regular.q@bs-dsvr05,regular.q@bs-dsvr08"
                        #queue="-q mpi04-ht.q" # bewi
                        echo "echo cd $dir_ > $fn_out; cd $dir_; echo $fn_bam >> $fn_out; hostname >> $fn_out; source `find ../ -name quality_statistics.sh`; source `find ../ -name paths.sh` ; source `find ../ -name utils.sh` ; quality_statistics $fn_bam $ext $stats_dir $gtf_file $out_dir $chromInfo $tsv_dir " | qsub $queue -o $fn_out -e $fn_err -N quality$cnt -cwd


                done
                exit 0
        fi

        echo run quality_statistics

        mkdir -p $out_dir

        for fn_bam in `echo $bamFileS`; do

                echo current bam $fn_bam
		sample_name=`basename ${fn_bam%$ext}`
                echo "sample_name = $sample_name"

		quality_statistics=$out_dir/${sample_name}_read_statistics.txt
		if [ ! -f $quality_statistics ]; then
			
			touch $quality_statistics
			# find which original bam files were used to be merged to this bam file
			orig_tsv_files=""
			for tsv in `find $tsv_dir -name *.tsv`
			do
				ext_sample_name=`extract $tsv EXTERNAL_SAMPLE_NAME`
				if [ "$ext_sample_name" == "$sample_name" ]; then
					orig_tsv_files="$orig_tsv_files  `basename ${tsv%.tsv} | sed 's/_UnP[12]//'`"
				fi
			done
			echo $orig_tsv_files

	
			local summary=$out_dir/`basename ${fn_bam%.bam}`_fastqc/summary.txt
	                if [ ! -f $summary ]; then
	                        echo "$fastqc $fn_bam -o $out_dir -f bam"
	                        $fastqc $fn_bam -o $out_dir -f bam
	                fi
	                cat $summary

			total_nr=0
			total_nr_cliptrim=0
			total_nr_reads_fastq=0
			total_nr_reads_fastq_cliptrim=0			

			for tsv in `echo $orig_tsv_files | tr ' ' '\n' | sort | uniq`; do # it happened that the same tsv file is listed several times
			
				stats_file=$stats_dir/`basename ${tsv%.tsv} | sed 's/_UnP[12]//'`.stats
				echo "current stats file = $stats_file"
	
				nr_reads_fastq_R1=$(cat $stats_file | awk 'NR==2')
				nr_reads_fastq_R2=$(cat $stats_file | awk 'NR==4')
				echo "# reads in R1 raw = $nr_reads_fastq_R1"
				echo "# reads in R2 raw = $nr_reads_fastq_R2"

				total_nr_reads_fastq=`echo "( $nr_reads_fastq_R1 + $nr_reads_fastq_R2 )" | bc`
				echo "This makes together = $total_nr_reads_fastq"

				total_nr=$((total_nr + total_nr_reads_fastq))
				echo "Therefore, the total # of raw reads from the lanes is now at = $total_nr"				

				nr_reads_fastq_cliptrim_R1=$(cat $stats_file | awk 'NR==7')
				echo "# reads after cliptrim in R1 = $nr_reads_fastq_cliptrim_R1"
				nr_reads_fastq_cliptrim_R2=$(cat $stats_file | awk 'NR==9')
				echo "# reads after cliptrim in R2 = $nr_reads_fastq_cliptrim_R2"
				nr_reads_fastq_cliptrim_R1_UnP=$(cat $stats_file | awk 'NR==11')
				echo "# reads after cliptrim in R1 UnP = $nr_reads_fastq_cliptrim_R1_UnP"
				nr_reads_fastq_cliptrim_R2_UnP=$(cat $stats_file | awk 'NR==13')
				echo "# reads after cliptrim in R2 UnP = $nr_reads_fastq_cliptrim_R2_UnP"

				total_nr_reads_fastq_cliptrim=`echo "($nr_reads_fastq_cliptrim_R1 + $nr_reads_fastq_cliptrim_R2 + $nr_reads_fastq_cliptrim_R1_UnP + $nr_reads_fastq_cliptrim_R2_UnP )" | bc`
				echo "Therefore, the total number of reads after cliptrim is  = $total_nr_reads_fastq_cliptrim"
				total_nr_cliptrim=$((total_nr_cliptrim+total_nr_reads_fastq_cliptrim))			
				echo "The total number of reads after cliptrim from all previous lanes and this one combined is = $total_nr_cliptrim"	

			done

			# count reads in bam
			echo "counting reads in bam file $fn_bam:"
			cnt_bam=`$samtools view -F 4 $fn_bam | wc -l` # only those that are really aligned
			echo $cnt_bam

			NOTrmdup_bam=$(echo $fn_bam | sed 's/.rmdup//')
			if [ -f $NOTrmdup_bam ]; then
				echo "counting reads in bam file $NOTrmdup_bam:"	
				cnt_NOTrmdup_bam=`$samtools view -F 4 $NOTrmdup_bam | wc -l `
				echo $cnt_NOTrmdup_bam
				cnt_NOTrmdup_bam_percentage=`echo "scale=3; ($cnt_NOTrmdup_bam/$total_nr_cliptrim)*100" | bc`
			fi

			# percentages	
			clip_trim_percentage=`echo "scale=3; ($total_nr_cliptrim/$total_nr)*100 " | bc`
			cnt_bam_percentage=`echo "scale =3; ($cnt_bam/$total_nr_cliptrim) * 100" | bc`

			echo -e "$total_nr\ttotal number of reads in the fastq files" > $quality_statistics
			echo -e "$total_nr_cliptrim\ttotal number of reads in the clipped and trimmed fastq files" >> $quality_statistics
			echo -e "$clip_trim_percentage\tpercentage of reads that survived the clipping and trimming" >> $quality_statistics
			if [ -f $NOTrmdup_bam ]; then
				echo -e "$cnt_NOTrmdup_bam\ttotal number of reads aligned (including PCR duplicates)" >> $quality_statistics
				echo -e "$cnt_NOTrmdup_bam_percentage\tpercentage of reads aligned (including PCR duplicates)" >> $quality_statistics
			fi
			echo -e "$cnt_bam\ttotal number of reads aligned" >> $quality_statistics
			echo -e "$cnt_bam_percentage\tpercentage of reads aligned" >> $quality_statistics

			echo -e "$total_nr_reads_fastq\ttotal number of reads in the fastq file"
			echo -e "$total_nr_reads_fastq_cliptrim\ttotal number of reads in the clipped and trimmed fastq files"
			echo -e "$clip_trim_percentage\tpercentage of reads that survived the clipping and trimming"
			if [ -f $NOTrmdup_bam ]; then
                                echo -e "$cnt_NOTrmdup_bam\ttotal number of reads aligned (including PCR duplicates)"
                                echo -e "$cnt_NOTrmdup_bam_percentage\tpercentage of reads aligned (including PCR duplicates)"
                        fi
			echo -e "$cnt_bam\ttotal number of reads aligned"
			echo -e "$cnt_bam_percentage\tpercentage of reads aligned"

		fi
		
		# coverage
		mkdir -p $out_dir/coverage
		TempPerBaseCoverage=/tmp/perBaseCoverage_${sample_name}.txt
		PerBaseCoverage=$out_dir/coverage/perBaseCoverage_${sample_name}.txt.gz
		PerBaseCoverageR=$out_dir/coverage/perBaseCoverage_${sample_name}_R.txt
		TempNOIntersect=/tmp/NOintersect_${sample_name}.bam
		bedGraph=$out_dir/coverage/${sample_name}.bedGraph
		bigWig=$out_dir/coverage/${sample_name}.bw
		

		echo "These are the files that are being generated for the coverage analysis:"
		echo "TempPerBaseCoverage = $TempPerBaseCoverage"
		echo "PerBaseCoverage = $PerBaseCoverage"
		echo "PerBaseCoverageR = $PerBaseCoverageR"
		echo "TempNOIntersect = $TempNOIntersect"
		echo "bedGraph = $bedGraph"
		echo "bigWig = $bigWig"

		if [ ! -f $PerBaseCoverageR -a ! -f $PerBaseCoverage ]; then
			#__ For each base in the gene track, compute the coverage
			echo "$samtools view -b -F 4 $fn_bam | $coverageBed -abam stdin -b $gtf_file -d > $TempPerBaseCoverage"
			$samtools view -b -F 4 $fn_bam | $coverageBed -abam stdin -b $gtf_file -d > $TempPerBaseCoverage
			echo "cat $TempPerBaseCoverage | gzip -9 > $PerBaseCoverage"
			cat $TempPerBaseCoverage | gzip -9 > $PerBaseCoverage
			echo "zcat $PerBaseCoverage | awk '{print $NF}' > $PerBaseCoverageR"
			zcat $PerBaseCoverage | awk '{print $NF}' > $PerBaseCoverageR
			rm $TempPerBaseCoverage
		
			`find ../ -name coverageR_quality_stats.R` $PerBaseCoverageR 
			`find ../ -name coverageR_quality_stats.R` $PerBaseCoverageR >> $quality_statistics
		fi

		if [ ! -f $TempNOIntersect ]; then
			echo "$samtools view -b -F 4 $fn_bam | $intersectBed -abam stdin -b $gtf_file -v > $TempNOIntersect"
			$samtools view -b -F 4 $fn_bam | $intersectBed -abam stdin -b $gtf_file -v > $TempNOIntersect  # alignments in this bam file do not overlap the gene tracks
        		cnt_bam_NO=`$samtools view -F 4 $TempNOIntersect | wc -l` # number of alignments NOT overlapping the gene tracks
		        cnt_bam_YES=$(( $cnt_bam - $cnt_bam_NO )) # number of alignments overlapping the gene tracks
		        echo "For ${sample_name}: $cnt_bam alignments in total; $cnt_bam_YES of them overlap the gene tracks; $cnt_bam_NO do not overlap gene tracks"
		        #__ percentage
		        cnt_percentage_YES=`echo "scale =3; ($cnt_bam_YES/$cnt_bam) * 100" | bc`
		        echo -e "$cnt_percentage_YES\tpercentage of the alignments that overlap the gene tracks" >> $quality_statistics
			echo -e "$cnt_percentage_YES\tpercentage of the alignments that overlap the gene tracks"
		        rm $TempNOIntersect
		fi

		if [ ! -f $bedGraph ]; then
			echo "$genomeCoverageBed  -ibam $fn_bam -bg -split -g $chromInfo > $bedGraph"
			$genomeCoverageBed  -ibam $fn_bam -bg -split -g $chromInfo > $bedGraph
		fi
		if [ ! -f $bigWig ]; then
			echo "$wigToBigWig $bedGraph $chromInfo $bigWig"
			$wigToBigWig $bedGraph $chromInfo $bigWig
		fi
	done
}


