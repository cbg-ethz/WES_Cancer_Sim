#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`
source `find ${gitDir} -name paths.sh`
source `find ${gitDir} -name submit.sh`
source `find ${gitDir} -name genome.sh`
source `find ${gitDir} -name bowtie2.sh`
source `find ${gitDir} -name tools.sh`
source `find ${gitDir} -name utils.sh`
ref=$fn_genome
run_wessim=`find ${gitDir} -name run_wessim.sh`

SGE_stuff="export SGE_ROOT=/usr/local/grid/soge; export BSSE_SGE=1; export SGE_CELL=bsse; export SGE_ARCH=lx-amd64; export SGE_CLUSTER_NAME=bsse-SoGE; export PATH=/usr/local/grid/soge/bin/lx-amd64/:\$PATH"

design=S04380219 #:v5 +UTR Nov 2012
fn_probes=${pathToBed}/$design/${design}_Probes.txt

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

if false ; then
for percentage in 0 0.25 0.33 0.5
do

    fasta_dir=${fasta_sim_dir}/final_vcfs_${percentage}
    out_dir=$fasta_dir/tree_blat_i${blat_min_identity}_s${blat_min_score}
    mkdir -p $out_dir
    fn_abundance=$out_dir/abundance.txt 
    
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
        
        echo $vec > $out_dir/abundance.txt 
        echo $cnt >> $out_dir/abundance.txt 
    fi

    echo $vec 
    port=6666
    cnt=7
    for num in $vec _NO_final; do
        for g in 0 1; do
            if [ $num == _NO_final ]; then 
                tag=${g}_NO_final
                fn_fasta=$out_dir/$tag.fa
                command="${run_wessim} $fn_fasta $fn_probes $out_dir/$tag 50000000 --num-cpus 6 --blat-port $port --blat-min-identity $blat_min_identity --blat-min-score $blat_min_score"
            else
                tag=${cnt}_${g}_TU_final
                fn_fasta=$out_dir/$tag.fa
                command="${run_wessim} $fn_fasta $fn_probes $out_dir/$tag $num --num-cpus 4 --blat-port $port --blat-min-identity $blat_min_identity --blat-min-score $blat_min_score"
            fi

            if [ ! -f "$fn_fasta" ]; then 
                echo "ln -s $fasta_dir/$tag.fa $fn_fasta"
                ln -s $fasta_dir/$tag.fa $fn_fasta
            fi

            port=$((port+1))

            #if true ; then	
            if false ; then	
                echo run locally
                $command
            else
                mkdir -p $out_dir/$tag
                fn_out=$out_dir/$tag/log.o
                if [ ! -f $fn_out ]; then 
                    echo "submit \"source ~/.bashrc; $command\" --tag CS_$tag --fn-out $fn_out"
                    submit "source ~/.bashrc; $command" --tag CS_$tag --fn-out $fn_out
                fi
            fi
        done
        cnt=$((cnt +1))
    done
done
wait_for_jobs CS_
fi

for percentage in 0.25 0 0.33 0.5
do

    fasta_dir=${fasta_sim_dir}/final_vcfs_${percentage}
    out_dir=$fasta_dir/tree_blat_i${blat_min_identity}_s${blat_min_score}
    mkdir -p $out_dir
    fn_abundance=$out_dir/abundance.txt 
    
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
        
        echo $vec > $out_dir/abundance_test.txt 
        echo $cnt >> $out_dir/abundance_test.txt 
    fi

    echo $vec 

    fasta_dir=${fasta_sim_dir}/final_vcfs_${percentage} 
    out_dir=$fasta_dir/tree_blat_i${blat_min_identity}_s${blat_min_score}

    if true; then 
        for bowtie_opts in "--sensitivity --sensitive --k 1 " "--sensitivity --very-sensitive --k 1" "--sensitivity --very-sensitive --k 100" "--sensitivity --very-sensitive --k 20 " ; do

            my_cnt=7
            sn=vsn
            if [ -z "`echo $bowtie_opts | grep very-sensitive`" ]; then 
                sn=sn
            fi

            align_dir=$out_dir/alignments_${sn}_k`echo $bowtie_opts | cut -f 2 -d "k" | cut -f 2 -d " "`

            #for num in $vec _NO_final; do
            for num in $vec; do
                for g in 0 1; do
                    if [ $num == _NO_final ]; then 
                        tag=${g}_NO_final
                    else
                        tag=${my_cnt}_${g}_TU_final
                    fi
                    fastq1=$out_dir/$tag/result_1.fastq.gz
                    fastq2=$out_dir/$tag/result_2.fastq.gz
                    out_prefix=$tag.tsv
                    command=" source `find ${dir_} -name bowtie2.sh`; source `find ${dir_} -name genome.sh`; source `find ${dir_} -name tools.sh`; source `find ${dir_} -name utils.sh`; run_botie2 $ref PAIRED_END $fastq1 $fastq2 $align_dir $out_prefix  --sample $tag $bowtie_opts"

                    if false ; then	
                        echo run locally
                        $command
                    else
                        #mkdir -p $out_dir/$tag
                        #fn_out=$out_dir/$tag/log_align.o
                        fn_out=$align_dir/logs/align_$tag.o
                        if [ ! -f $fn_out ]; then 
                            echo "submit \"$command\" --tag align_$tag --fn-out $fn_out"
                            submit "$command" --tag align_$tag --fn-out $fn_out
                        fi
                    fi
                done
                my_cnt=$((my_cnt +1))
            done
        done
    fi
done
wait_for_jobs align_
echo alignments done

# run variant caller script
################################################################################
for align_dir in $out_dir/alignments*; do
	echo ${gitDir}/sim_cancer/simulate_WES_cancerSample/process_alignments.sh $align_dir/bam/ 
	${gitDir}/sim_cancer/simulate_WES_cancerSample/process_alignments.sh $align_dir/bam/ & 
	#echo ./process_true_alignments.sh $align_dir/bam/
	${gitDir}/sim_cancer/simulate_WES_cancerSample/process_true_alignments.sh $align_dir/bam/ & 
	echo ${gitDir}/sim_cancer/var_callers/run_variant_callers.sh $align_dir/bam/ $align_dir/variant_callers/  
	${gitDir}/sim_cancer/var_callers/run_variant_callers.sh $align_dir/bam/ $align_dir/variant_callers/ &
done
exit 0

if false ; then 
	for fn_bam in $align_dir/bam/*sorted.bam
	do
		fn_best=${fn_bam%.bam}.best.bam
		echo "$samtools view -b -F 256 $fn_bam > $fn_best"
		$samtools view -b -F 256 $fn_bam > $fn_best
	done
	exit 0
fi


