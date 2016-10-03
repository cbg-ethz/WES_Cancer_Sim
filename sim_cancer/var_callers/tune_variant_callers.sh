#/bin/bash -l

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`

#cd $SNV_DIR
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


# parse arguments
outputAlignDir=$1
outputDir=$2
bam_ending=$3
if [ -z "$3" ]; then
	bam_ending=perc.bam
fi


################################################################################################################################
#
## Run pipeline
#

### set the tag: this will determine the name-tag of all the output folders
#tag=default
tag=tuned

##queue
queue=mpi04-ht.q

# set whether or not local realignment around indels should be done (in this mode, only cont=20 and perc=50 will be analyzed)
#local_realign=true
local_realign=false

bamDir=$outputAlignDir
variantDir=$outputDir/variants_${tag}
#qualityOut=$outputDir/quality_control_${tag}
mkdir -p $bamDir
mkdir -p $variantDir
#mkdir -p $qualityOut

# process the alingments and generate the various coverage levels/ contamination files
num_bams=`ls -lh $bamDir/TU.wCont[0-9]0.final.RG.*$bam_ending $bamDir/NO_final.RG.*$bam_ending | wc -l`
if [ $num_bams -lt 13 ]; then
	echo "start process alignments in 60 seconds"
	sleep 60
	command=`find $currScriptDir/../ -name process_alignments.sh`
	eval $command
fi
# cleanup intermediate bam files
#rm $bamDir/*.sorted.woMM.bam
#rm $bamDir/*.sorted.woMM.karyotypic.bam
#rm $bamDir/NO_final.[0-9]*.bam*
#rm $bamDir/TU.wCont[0-9]0.final.bam
#rm $bamDir/NO_final.bam $bamDir/TU_final.bam


### perform quality control 
##log_dir=$qualityOut/logs/
##for curr_bam in $bamDir/TU.wCont[0-9]0.final.RG.*$bam_ending $bamDir/NO_final.RG.*$bam_ending
##do
##	echo $curr_bam
##	fn_out=${log_dir}cov_job_`basename ${curr_bam%.bam}`.o
##
##	if [ -f $fn_out ]; then
##		continue;
##	fi
##	command="source `find $SNV_DIR -name coverage.sh` ; source `find $SNV_DIR -name paths.sh` ; coverage $curr_bam $bed_file $qualityOut"
##	submit "$command" --tag cov_stats --log-dir $log_dir
##done

####################################
#### call variants # the different "percent files"


# set the percentage (=coverage) and the contamination level
perc=50
cont=20

if $local_realign; then
	echo "perform local realignment for the bam files with normal contamination $cont and percentage $perc"
	for TU_bam in `find $bamDir -name TU.wCont[0-9]0.final.RG.*$bam_ending | grep ${perc}perc | grep wCont${cont}`
	do
		NO_bam=`echo $TU_bam | sed 's/TU/NO/' | sed 's/.wCont[0-9]0./_/'`;
		if [ ! -f $NO_bam ]; then
			NO_bam=$(echo $NO_bam | sed 's/final.RG./final.RG.true./') # the normal true bams have an extra "true" in their name
		fi

		log_dir=$bamDir/logs/
		fn_out=$log_dir/loc_realign_`basename ${TU_bam%.bam}`.o
		if [ -f $fn_out ]; then
			echo "$fn_out already exists."
		else
			outBam=${TU_bam%.bam}.locRealign.bam
			command="source `find $SNV_DIR -name loc_realign.sh` ; source `find $SNV_DIR -name paths.sh`; loc_realign $TU_bam $outBam $fn_genome --create-index --check-fexist "
			submit "$command" --fn-out $fn_out --tag locRealign`basename $TU_bam` --work-dir $SNV_DIR --queue $queue
		fi

		fn_out=$log_dir/loc_realign_`basename ${NO_bam%.bam}`.o
		if [ -f $fn_out ]; then
			echo "$fn_out already exists."
		else
			outBamNO=${NO_bam%.bam}.locRealign.bam
			command="source `find $SNV_DIR -name loc_realign.sh`; source `find $SNV_DIR -name paths.sh`; loc_realign $NO_bam $outBamNO $fn_genome --create-index --check-fexist "
			submit "$command" --fn-out $fn_out --tag locRealign`basename $NO_bam` --work-dir $SNV_DIR --queue $queue
		fi
	done
	wait_for_jobs locRealign
fi
 

echo start tuning the variant callers... 
for TU_bam in `find $bamDir -name TU.wCont[0-9]0.final.RG.*$bam_ending | grep ${perc}perc | grep wCont${cont}`
do
	NO_bam=`echo $TU_bam | sed 's/TU/NO/' | sed 's/.wCont[0-9]0./_/'`;
	if [ ! -f $NO_bam ]; then
		NO_bam=$(echo $NO_bam | sed 's/final.RG./final.RG.true./') # the normal true bams have an extra "true" in their name
	fi
	echo "TU_bam = $TU_bam"
	echo "NO_bam = $NO_bam"

	if $local_realign; then
		TU_bam=${TU_bam%.bam}.locRealign.bam
		NO_bam=${NO_bam%.bam}.locRealign.bam
	fi

	sample_name=TU.${perc}perc
	log_dir=$variantDir/logs/
	mkdir -p $log_dir

	
	#samtools_1_2
	################################################################################
	declare -a sam_opts=('-C 20' '-C 30' '-C 40' '-C 50' '-C 70' '-B' '-d 500' '-Q 0' '-Q 25' '-q 5' '-q 20' '-C 100' '-C 150' '-C 200' '-C 250' '-C 500' 'withE');
	for opt_cnt in `seq 1 ${#sam_opts[@]}`
	do
		echo current option: ${sam_opts[opt_cnt-1]}
		fn_out=$log_dir/samvar_1_2_`basename ${TU_bam%.bam}`option_`echo ${sam_opts[opt_cnt-1]} | sed 's/-//g' | sed 's/ //g'`.o
		if [ -f $fn_out ]; then
				echo "$fn_out already exists."
		else
			if [ "${sam_opts[opt_cnt-1]}" == "withE" ]; then
				command="source `find $SNV_DIR -name samtools_1_2_variants.sh` ; source `find $SNV_DIR -name paths.sh`; samtools_1_2_variants $TU_bam xxx $fn_genome $variantDir $NO_bam "
				submit "$command" --fn-out $fn_out --tag sam_1_2`basename ${TU_bam}` --work-dir $SNV_DIR --queue $queue
			else
				command="source `find $SNV_DIR -name samtools_1_2_variants.sh` ; source `find $SNV_DIR -name paths.sh`; samtools_1_2_variants $TU_bam xxx $fn_genome $variantDir $NO_bam --opts '${sam_opts[opt_cnt-1]} --no-E' "
				submit "$command" --fn-out $fn_out --tag sam_1_2`basename ${TU_bam}` --work-dir $SNV_DIR --queue $queue
			fi
		fi

		# also for the paired calls TU+NO together:
		fn_out=$log_dir/samvar_1_2_paired_TU_NO_`basename ${TU_bam%.bam}`option_`echo ${sam_opts[opt_cnt-1]} | sed 's/-//g' | sed 's/ //g'`.o
		if [ -f $fn_out ]; then
				echo "$fn_out already exists."
		else
			if [ "${sam_opts[opt_cnt-1]}" == "withE" ]; then
				command="source `find $SNV_DIR -name samtools_1_2_variants.sh` ; source `find $SNV_DIR -name paths.sh`; samtools_1_2_variants $TU_bam xxx $fn_genome $variantDir $NO_bam --paired"
				submit "$command" --fn-out $fn_out --tag sam_1_2`basename ${TU_bam}` --work-dir $SNV_DIR --queue $queue
			else
				command="source `find $SNV_DIR -name samtools_1_2_variants.sh` ; source `find $SNV_DIR -name paths.sh`; samtools_1_2_variants $TU_bam xxx $fn_genome $variantDir $NO_bam --opts '${sam_opts[opt_cnt-1]} --no-E' --paired "
				submit "$command" --fn-out $fn_out --tag sam_1_2`basename ${TU_bam}` --work-dir $SNV_DIR --queue $queue
			fi
		fi
	done

	#deepSNV
	################################################################################
	declare -a deep_opts=('--minBaseQ 0' '--minBaseQ 13' '--overdispersion 1000' '--overdispersion 500' '--overdispersion 200' )
	for opt_cnt in `seq 1 ${#deep_opts[@]}`
	do
		echo current option is: ${deep_opts[opt_cnt-1]}
		fn_out=${log_dir}deepSNV_job_`basename ${TU_bam%.bam}`option_`echo ${deep_opts[opt_cnt-1]} | sed 's/-//g' | sed 's/ //g'`.o
		if [ -f $fn_out ]; then
			echo "$fn_out already exists."
		else
			command="source `find $SNV_DIR -name deepSNV_automated_variants.sh` ; source `find $SNV_DIR -name paths.sh` ; deepSNV_automated_variants $bed_file $TU_bam $NO_bam $variantDir --also-vcf true ${deep_opts[opt_cnt-1]} "
			submit "$command" --fn-out $fn_out --tag deepSNV`basename ${TU_bam}` --work-dir $SNV_DIR --queue $queue
		fi
	done


	#jointSNVMix2
	################################################################################
	declare -a jsm_opts1=('0' '13' '25')
	declare -a jsm_opts2=('0' '5' '20')
	for opt1_cnt in `seq 1 ${#jsm_opts1[@]}`
	do
		for opt2_cnt in `seq 1 ${#jsm_opts2[@]}`
		do
			#echo ${jsm_opts1[opt1_cnt-1]}
			#echo ${jsm_opts2[opt2_cnt-1]}
			if [ ${jsm_opts1[opt1_cnt-1]} == 0 -a ${jsm_opts2[opt2_cnt-1]} == 0 ]; then
				continue
			fi
			curr_opt1=" --min_base_qual ${jsm_opts1[opt1_cnt-1]} "
			curr_opt2=" --min_map_qual ${jsm_opts2[opt2_cnt-1]} "

			fn_out=${log_dir}jointSNVMix2_job_`basename ${TU_bam%.bam}`_option_`echo $curr_opt1 | sed 's/-//g' | sed 's/ //g'`_`echo $curr_opt2 | sed 's/-//g' | sed 's/ //g'`.o
			if [ -f $fn_out ]; then
				echo "$fn_out already exists."
			else
				command="source `find $SNV_DIR -name jointSNVMix2_variants.sh` ; source `find $SNV_DIR -name paths.sh` ; jointSNVMix2_variants $TU_bam '' $fn_genome $variantDir $NO_bam $bed_file $curr_opt1 $curr_opt2 "
				submit "$command" --fn-out $fn_out --tag jointSNVMix2$perc --work-dir $SNV_DIR --queue $queue
				echo $command
			fi
		done
	done

	#varscan2
	################################################################################ 
	declare -a varscan_opts=('--min-var-freq 0.0001' '--min-var-freq 0.005' '--min-var-freq 0.01'  '--min-var-freq 0.02' '--min-var-freq 0.05' '--tumor-purity 0.80' '--strand-filter 1' '--min-tumor-freq 0.02' '--min-tumor-freq 0.05') # the last two are for processSomatic
	for opt_cnt in `seq 1 ${#varscan_opts[@]}`
	do
		curr_opt=${varscan_opts[opt_cnt-1]}
		fn_out=${log_dir}varscan2_job_`basename ${TU_bam%.bam}`_option_`echo $curr_opt | sed 's/-//g' | sed 's/ //g'`.o
		if [ -f $fn_out ]; then
			echo "$fn_out already exists."
		else
			command="source `find $SNV_DIR -name varscan2_variants.sh` ; source `find $SNV_DIR -name paths.sh` ; varscan2_variants $TU_bam xxx $fn_genome $variantDir $NO_bam $curr_opt "
			submit "$command" --fn-out $fn_out --tag varscan2$perc --work-dir $SNV_DIR --queue $queue
			echo $command
		fi
	done

	# sinvict
	################################################################################
	declare -a sinvict_opts=('--qscore-cutoff 20' '--qscore-cutoff 40' '--qscore-cutoff 60' '--qscore-cutoff 80' )
	for opt_cnt in `seq 1 ${#sinvict_opts[@]}`
	do
		curr_opt=${sinvict_opts[opt_cnt-1]}
		fn_out=${log_dir}sinvict_job_`basename ${TU_bam%.bam}`_option_`echo $curr_opt | sed 's/-//g' | sed 's/ //g'`.o
		if [ -f $fn_out ]; then
			echo "$fn_out already exists."
		else
			command="source `find $SNV_DIR -name sinvict_variants.sh`; source `find $SNV_DIR -name paths.sh` ; sinvict_variants $TU_bam $fn_genome $variantDir $NO_bam $bed_file $curr_opt "
			submit "$command" --fn-out $fn_out --tag sinvict$perc --work-dir $SNV_DIR --queue $queue
		fi
	done

done
exit 0




