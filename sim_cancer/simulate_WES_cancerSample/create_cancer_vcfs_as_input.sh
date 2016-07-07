#/bin/bash -l

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`
source `find $SNV_DIR -name paths.sh`
source `find $SNV_DIR -name utils.sh`
source `find $SNV_DIR -name submit.sh`
source `find $SNV_DIR -name merge_samples_tsv.sh`
source `find $SNV_DIR -name fix_bam_header.sh`
source `find $SNV_DIR -name quality_statistics.sh`
source `find $SNV_DIR -name picard_rmdup.sh`
source `find $SNV_DIR -name freebayes_variants.sh`
source `find $SNV_DIR -name varscan2_variants.sh`
source `find $SNV_DIR -name muTect_variants.sh`
source `find $SNV_DIR -name samtools_variants.sh`
source `find $SNV_DIR -name samtools_1_2_variants.sh`
source `find $SNV_DIR -name coverage.sh`
source `find $SNV_DIR -name gatkUnifiedGT_variants.sh`
source `find $SNV_DIR -name somaticSniper_variants.sh`
source `find $SNV_DIR -name samtools_rmMM.sh`
source `find $SNV_DIR -name sort_bam.sh`
source `find $SNV_DIR -name subsetSam.sh`
source `find $SNV_DIR -name replaceReadGroups.sh`
source `find $SNV_DIR -name merge_bams.sh`
source `find $SNV_DIR -name loc_realign.sh`
source `find $SNV_DIR -name gatkHPCaller_variants.sh`

bamDir=$1
if [ -z "$1" ]; then
	echo "Please provide directory with bams"
	exit 0
fi

fn_genome=$humanRefKaryo

log_dir=$bamDir/logs
mkdir -p $log_dir
outDir=$bamDir/freebayes_out
mkdir -p $outDir
declare -a freebayes_opts1=('2' '3' '4' '5' '6')
declare -a freebayes_opts2=('0.01' '0.05' '0.1' '0.2')
pooled_continuous="pooled_continuous"
min_mapping_quality=0

# default is --min-alternate-count 2 --min-alternate-fraction 0.2

for fn_bam in `find $bamDir -name *.bam`
do
	for opts1_cnt in `seq 1 ${#freebayes_opts1[@]}`
	do
		for opts_2_cnt in `seq 1 ${#freebayes_opts2[@]}`
		do
			if [ ${freebayes_opts1[opts1_cnt-1]} == 2 -a `echo ${freebayes_opts2[opts_2_cnt-1]}'=='0.01 | bc -l` == 1 ]; then
				continue
			fi
			curr_opt1=" --min-alternate-count ${freebayes_opts1[opts1_cnt-1]} "
			curr_opt2=" --min-alternate-fraction ${freebayes_opts2[opts_2_cnt-1]} "

			fn_out=$log_dir/freebayes_`basename ${fn_bam%.bam}`_option_`echo $curr_opt1 | sed 's/-//g' | sed 's/ //g'`_`echo $curr_opt2 | sed 's/-//g' | sed 's/ //g'`.o
			if [ -f $fn_out ]; then
				echo "$fn_out already exists."
			else
				command="source `find $SNV_DIR -name freebayes_variants.sh` ; source `find $SNV_DIR -name paths.sh`; freebayes_variants $fn_bam '' $fn_genome $outDir $curr_opt1 $curr_opt2 --pooled-continuous $pooled_continuous --min-mapping-quality $min_mapping_quality "
				echo $command
				submit "$command" --fn-out $fn_out --tag freebayes_${freebayes_opts1[opts1_cnt-1]}_${freebayes_opts2[opts_2_cnt-1]}_`basename $fn_bam` --work-dir $SNV_DIR --queue mpi04-ht.q
			fi
		done
	done
done

