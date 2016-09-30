#!/bin/bash

# Description:
# This script computes the coverage of a bam file. Output file contains: minimum, maximum, mean, and median coverage of the genes/bed file;
# % of bases in the genes/bed file that are not covered at all, covered with at least 20 reads;
# % of reads that overlap the genes/bed file. It also executes fastqc.

# Inputs:
# BamFiles = One or several bam files for which the quality control should be conducted.
# gtf_file_genes = the gtf file that contains only the gene boundaries (exons, introns, etc. removed) OR the BED file with the targeted regions of the exome kit
# output_dir = the output folder for all the quality control files

# Outputs:
# Within the output folder there will be the following directories:
# coverage = for each sample the files used to determine the coverage; most of them generated with BEDtools
# <bam_file%.bam>_fastqc = the folder which serves as output from the program fastqc
# In the output folder there will be the following files:
# <bam_file%.bam>_cov_stats.txt = containing a summary of all the coverage information


function coverage()
{

        if [ -z "$1" -o -z "$2" -o -z "$3" ]; then
                echo "usage coverage <BamFiles> <gtf_file_genes> <output_dir>"
                exit -1;
        fi

        local bamFileS=$1; shift;
        local gtf_file=$1; shift;
        local out_dir=$1; shift;

        echo "Running coverage using the following inputs:"
        echo "bamFileS = $bamFileS"
	echo "gtf_file_genes = $gtf_file"
        echo "out_dir = $out_dir"

        echo run coverage

        mkdir -p $out_dir

        for fn_bam in `echo $bamFileS`; do

                echo current bam $fn_bam

		coverage_stats=$out_dir/`basename ${fn_bam%.bam}`_cov_stats.txt
		if [ ! -f $coverage_stats ]; then
			
			touch $coverage_stats
	
			local summary=$out_dir/`basename ${fn_bam%.bam}`_fastqc/summary.txt
	                if [ ! -f $summary ]; then
	                        echo "$fastqc $fn_bam -o $out_dir -f bam"
	                        $fastqc $fn_bam -o $out_dir -f bam
	                fi
	                cat $summary
		
			mkdir -p $out_dir/coverage
			TempPerBaseCoverage=/tmp/perBaseCoverage_`basename ${fn_bam%.bam}`.txt
			PerBaseCoverage=$out_dir/coverage/perBaseCoverage_`basename ${fn_bam%.bam}`.txt.gz
			PerBaseCoverageR=$out_dir/coverage/perBaseCoverage_`basename ${fn_bam%.bam}`_R.txt
			TempNOIntersect=/tmp/NOintersect_`basename ${fn_bam%.bam}`.bam
	
			echo "These are the files that are being generated for the coverage analysis:"
			echo "TempPerBaseCoverage = $TempPerBaseCoverage"
			echo "PerBaseCoverage = $PerBaseCoverage"
			echo "PerBaseCoverageR = $PerBaseCoverageR"
			echo "TempNOIntersect = $TempNOIntersect"
	
			if [ ! -f $PerBaseCoverageR -a ! -f $PerBaseCoverage ]; then
				#__ For each base in the gene track, compute the coverage
				echo "$samtools view -b -F 4 $fn_bam | $coverageBed -abam stdin -b $gtf_file -d > $TempPerBaseCoverage"
				$samtools view -b -F 4 $fn_bam | $coverageBed -abam stdin -b $gtf_file -d > $TempPerBaseCoverage
				echo "cat $TempPerBaseCoverage | gzip -9 > $PerBaseCoverage"
				cat $TempPerBaseCoverage | gzip -9 > $PerBaseCoverage
				echo "zcat $PerBaseCoverage | awk '{print \$NF}' > $PerBaseCoverageR"
				zcat $PerBaseCoverage | awk '{print $NF}' > $PerBaseCoverageR
				rm $TempPerBaseCoverage
			
				`find $gitDir -name coverageR_quality_stats.R` $PerBaseCoverageR 
				`find $gitDir -name coverageR_quality_stats.R` $PerBaseCoverageR >> $coverage_stats
			fi
	
			if [ -f $TempNOIntersect ]; then
				rm $TempNOIntersect
			fi
			if [ ! -f $TempNOIntersect ]; then
				cnt_bam=`$samtools view -F 4 $fn_bam | wc -l`
				echo -e "$cnt_bam\tnumber of reads aligned"
				echo -e "$cnt_bam\tnumber of reads aligned" >> $coverage_stats
				echo "$samtools view -b -F 4 $fn_bam | $intersectBed -abam stdin -b $gtf_file -v > $TempNOIntersect"
				$samtools view -b -F 4 $fn_bam | $intersectBed -abam stdin -b $gtf_file -v > $TempNOIntersect  # alignments in this bam file do not overlap the gene tracks
	        		cnt_bam_NO=`$samtools view -F 4 $TempNOIntersect | wc -l` # number of alignments NOT overlapping the gene tracks
			        cnt_bam_YES=$(( $cnt_bam - $cnt_bam_NO )) # number of alignments overlapping the gene tracks
			        echo "For ${fn_bam%.bam}: $cnt_bam alignments in total; $cnt_bam_YES of them overlap the gene tracks; $cnt_bam_NO do not overlap gene tracks"
			        #__ percentage
			        cnt_percentage_YES=`echo "scale =3; ($cnt_bam_YES/$cnt_bam) * 100" | bc`
			        echo -e "$cnt_percentage_YES\tpercentage of the alignments that overlap the gene tracks" >> $coverage_stats
				echo -e "$cnt_percentage_YES\tpercentage of the alignments that overlap the gene tracks"
			        rm $TempNOIntersect
			fi
		else 
			echo "$coverage_stats already exists. Nothing to be done."
		fi

	done
}


