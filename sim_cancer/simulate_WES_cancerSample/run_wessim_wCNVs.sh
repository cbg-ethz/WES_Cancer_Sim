#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`
source `find ${gitDir} -name paths.sh`
source `find ${gitDir} -name submit.sh`
source `find ${gitDir} -name genome.sh`
source `find ${gitDir} -name bowtie2.sh`
source `find ${gitDir} -name utils.sh`
echo "bowtie = $bowtie"

ref=$fn_genome
run_wessim=`find ${gitDir} -name run_wessim.sh`

SGE_stuff="export SGE_ROOT=/usr/local/grid/soge; export BSSE_SGE=1; export SGE_CELL=bsse; export SGE_ARCH=lx-amd64; export SGE_CLUSTER_NAME=bsse-SoGE; export PATH=/usr/local/grid/soge/bin/lx-amd64/:\$PATH"

design=S04380219 #:v5 +UTR Nov 2012
fn_probes=${pathToBed}/$design/${design}_Probes.txt

fasta_dir_cnvs=$1
outDir_cnvs=$2

if [ -z "$1" -o -z "$2" ]; then
	echo "please provide fasta_dir_cnvs, outDir_cnvs"
	exit 0
fi

function chop_count()
{
	res=""
	for x in $@; do
		#x=$1
		#shift
		if [ $x -eq 0 ]; then 
			x1=0
			x2=0
		else
			x1=`echo $(( $RANDOM$RANDOM % $x))`;
			x2=$((x-$x1));
		fi
		res="$res $x1 $x2 "
	done
	echo $res
}

blat_min_identity=70
blat_min_score=60

mkdir -p $outDir_cnvs
fn_abundance=$outDir_cnvs/abundance.txt 
    
if [ -f $fn_abundance ]; then 

        vec=`head -n 1 $fn_abundance`

else
	num_reads_total=50000000
      
        # generate abundance vector
        #vec=`chop_count 10000`
        vec=`chop_count $num_reads_total`
        vec=`chop_count $vec`
        vec=`chop_count $vec`
        
        cnt=""
        for x in $vec; do val=`echo "scale=3; $x/$num_reads_total*100" | bc`; cnt="$cnt $val"; done
        echo cnt $cnt
        
        echo $vec > $outDir_cnvs/abundance.txt 
        echo $cnt >> $outDir_cnvs/abundance.txt 
fi


echo $vec 
port=6666
cnt=7
for num in $vec _NO_final; do
        for g in 0 1; do
            if [ $num == _NO_final ]; then 
                tag=${g}_NO_final
                fn_fasta=$fasta_dir_cnvs/$tag.fa
                command="${run_wessim} $fn_fasta $fn_probes $outDir_cnvs/$tag 50000000 --num-cpus 6 --blat-port $port --blat-min-identity $blat_min_identity --blat-min-score $blat_min_score"
            else
                tag=${cnt}_${g}_TU
                fn_fasta=$fasta_dir_cnvs/$tag.fa
                command="${run_wessim} $fn_fasta $fn_probes $outDir_cnvs/$tag $num --num-cpus 4 --blat-port $port --blat-min-identity $blat_min_identity --blat-min-score $blat_min_score"

		# in case there is a gain:
		tagGain=${cnt}_${g}_TU_gain
		fn_fastaGain=$fasta_dir_cnvs/$tagGain.fa
		commandGain="${run_wessim} $fn_fastaGain $fn_probes $outDir_cnvs/$tagGain $num --num-cpus 4 --blat-port $port --blat-min-identity $blat_min_identity --blat-min-score $blat_min_score"
            fi

            if [ ! -f "$fn_fasta" ]; then 
                echo "Error: $fn_fasta does not exist"
                exit 0
            fi

            port=$((port+1))

            #if true ; then	
            if false ; then	
                mkdir -p $outDir_cnvs/$tag
                fn_out=$outDir_cnvs/$tag/log.o
                if [ ! -f $fn_out ]; then 
                	echo run locally
                	eval $command > $fn_out 2>&1 
		fi

		# if there is a gain:
                if [ -f $fn_fastaGain ]; then
                        mkdir -p $outDir_cnvs/$tagGain
                        fn_out=$outDir_cnvs/$tagGain/log.o
                        if [ ! -f $fn_out ]; then
				eval $commandGain > $fn_out 2>&1 
                        fi
                fi
            else
                mkdir -p $outDir_cnvs/$tag
                fn_out=$outDir_cnvs/$tag/log.o
                if [ ! -f $fn_out ]; then 
                    echo "submit \"source ~/.bashrc; $command\" --tag CS_$tag --fn-out $fn_out"
                    submit "source ~/.bashrc; export PYTHONPATH=''; echo ${PYTHONPATH}; $command" --tag CS_$tag --fn-out $fn_out
                fi

		# if there is a gain:
		if [ -f $fn_fastaGain ]; then
			mkdir -p $outDir_cnvs/$tagGain
			fn_out=$outDir_cnvs/$tagGain/log.o
			if [ ! -f $fn_out ]; then
				echo "submit \"source ~/.bashrc; $commandGain\" --tag CS_$tagGain --fn-out $fn_out"
				submit "source ~/.bashrc; $commandGain" --tag CS_$tagGain --fn-out $fn_out
			fi
		fi
            fi
        done
        cnt=$((cnt +1))
    done
wait_for_jobs CS_


for bowtie_opts in  "--sensitivity --very-sensitive --k 20 " ### "--sensitivity --sensitive --k 1 " "--sensitivity --very-sensitive --k 1" "--sensitivity --very-sensitive --k 100" "--sensitivity --very-sensitive --k 20 "
do
	my_cnt=7
	sn=vsn
	if [ -z "`echo $bowtie_opts | grep very-sensitive`" ]; then 
		sn=sn
	fi

	align_dir=$outDir_cnvs/alignments_${sn}_k`echo $bowtie_opts | cut -f 2 -d "k" | cut -f 2 -d " "`

	for my_cnt in 7 8 9 10 11 12 13 14 _NO_final; do
        	for g in 0 1; do
        		if [ $my_cnt == _NO_final ]; then 
        	        	tag=${g}_NO_final
			else
        	        	tag=${my_cnt}_${g}_TU
			fi
        	        fastq1=$outDir_cnvs/$tag/result_1.fastq.gz
        	        fastq2=$outDir_cnvs/$tag/result_2.fastq.gz
        	        out_prefix=$tag.tsv
        	        command=" source `find ${dir_} -name paths.sh`; source `find ${dir_} -name bowtie2.sh`; source `find ${dir_} -name genome.sh`; source `find ${dir_} -name utils.sh`; echo bowtie = $bowtie; run_botie2 $ref PAIRED_END $fastq1 $fastq2 $align_dir $out_prefix  --sample $tag $bowtie_opts"
        	        
			if false; then	
       				# check first that the fastq file is not empty:
				FILESIZE=$(stat -c%s "$fastq1")
				if [ $FILESIZE -le 35 ]; then 
					echo "$fastq1 might be empty!"
				else
		                	#mkdir -p $out_dir/$tag
		                	#fn_out=$out_dir/$tag/log_align.o
		                	fn_out=$align_dir/logs/align_$tag.o
		                        if [ ! -f $fn_out ]; then 
			         		echo run locally
		                		eval $command > $fn_out 2>&1
					fi
				fi
        	        else
				# check first that the fastq file is not empty:
				FILESIZE=$(stat -c%s "$fastq1")
				if [ $FILESIZE -le 35 ]; then 
					echo "$fastq1 might be empty!"
				else
		                	#mkdir -p $out_dir/$tag
		                	#fn_out=$out_dir/$tag/log_align.o
		                	fn_out=$align_dir/logs/align_$tag.o
		                        if [ ! -f $fn_out ]; then 
		                            echo "submit source `find ${gitDir} -name paths.sh`; \"$command\" --tag align_$tag --fn-out $fn_out"
		                            submit "source `find ${gitDir} -name paths.sh`; $command" --tag align_$tag --fn-out $fn_out
					fi
        	                fi
        	        fi

			# in case there is a gain:
			tagGain=${my_cnt}_${g}_TU_gain
			fastq1=$outDir_cnvs/$tagGain/result_1.fastq.gz
			fastq2=$outDir_cnvs/$tagGain/result_2.fastq.gz
			out_prefix=$tagGain.tsv
			command=" source `find ${dir_} -name paths.sh`; source `find ${dir_} -name bowtie2.sh`; source `find ${dir_} -name genome.sh`; source `find ${dir_} -name utils.sh`; run_botie2 $ref PAIRED_END $fastq1 $fastq2 $align_dir $out_prefix --sample $tagGain $bowtie_opts"

			# if there is a gain:
			if [ -f $fastq1 ]; then
				# check first that the fastq file is not empty:
				FILESIZE=$(stat -c%s "$fastq1")
				if [ $FILESIZE -le 35 ]; then 
					echo "$fastq1 might be empty!"
				else
					if false; then
						fn_out=$align_dir/logs/align_$tagGain.o
        	                        	if [ ! -f $fn_out ]; then
        	                        	        echo run locally
        	                        	        eval $command > $fn_out 2>&1
        	                        	fi
					else
	
						fn_out=$align_dir/logs/align_$tagGain.o
						if [ ! -f $fn_out ]; then
							echo "submit source `find ${gitDir} -name paths.sh`; \"$command\" --tag align_$tagGain --fn-out $fn_out"
							submit "source `find ${gitDir} -name paths.sh`; $command" --tag align_$tagGain --fn-out $fn_out
						fi
					fi
				fi
			fi
		done
	done
done

wait_for_jobs align_
echo alignments done


# process alignments and run variant caller script
################################################################################
for align_dir in $outDir_cnvs/alignments*; do
	command=" ${gitDir}/sim_cancer/simulate_WES_cancerSample/process_alignments_wCNVs.sh $align_dir/bam/ $fasta_dir_cnvs "
	echo $command
	##eval $command
	
	command="${gitDir}/sim_cancer/simulate_WES_cancerSample/process_true_alignments_wCNVs.sh $align_dir/bam/ $fasta_dir_cnvs "
	echo $command
	eval $command
done
exit 0


