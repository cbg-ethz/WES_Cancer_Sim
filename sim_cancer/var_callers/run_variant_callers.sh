#/bin/bash -l

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`

echo $currScriptDir 
echo $SNV_DIR

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
source `find $SNV_DIR -name gatkHPCaller_variants.sh`


# parse arguments

#outputAlignDir=$1
bamDir=$1
outputDir=$2
bam_ending=$3
if [ -z "$3" ]; then
	bam_ending=perc.bam
fi
echo $bam_ending


fn_genome=$humanRefKaryo

################################################################################################################################
#
## Run pipeline
#

### set the tag: this will determine the name-tag of all the output folders
tag=default
#tag=repeatedSubsampling

# set whether or not local realignment around indels should be done (in this mode, only cont=20 and perc=50 will be analyzed)
#local_realign=true
local_realign=false

#queue=mpi01.q
queue=mpi04-ht.q

variantDir=$outputDir/variants_${tag}
qualityOut=$outputDir/quality_control_${tag}
mkdir -p $bamDir
mkdir -p $variantDir
mkdir -p $qualityOut

statusDir=$outputDir/status


# process the alingments and generate the various converage levels/ contamination files
echo $bamDir
echo $variantDir

if [ ! -f $bamDir/TU.wCont20.final.RG.50perc.bam ]; then 
	echo "start process alignments in 20 seconds"
	sleep 20
	echo "source ${gitDir}/simulate_WES_cancerSample/process_alignments.sh"
	source ${gitDir}/simulate_WES_cancerSample/process_alignments.sh
fi
# cleanup intermediate bam files
#rm $bamDir/*.sorted.woMM.bam
#rm $bamDir/*.sorted.woMM.karyotypic.bam
#rm $bamDir/NO_final.[0-9]*.bam*
#rm $bamDir/TU.wCont[0-9]0.final.bam
#rm $bamDir/NO_final.bam $bamDir/TU_final.bam


### perform quality control 
if false; then
	log_dir=$qualityOut/logs/
	for curr_bam in $bamDir/TU.wCont[0-9]0.final.RG.*$bam_ending $bamDir/NO_final.RG.*$bam_ending
	do
		echo $curr_bam
		fn_out=${log_dir}cov_job_`basename ${curr_bam%.bam}`.o
	
		if [ -f $fn_out ]; then
			continue;
		fi
		command="source `find $SNV_DIR -name coverage.sh` ; source `find $SNV_DIR -name paths.sh` ; coverage $curr_bam $bed_file $qualityOut"
		submit "source `find $SNV_DIR -name paths.sh`; source `find $SNV_DIR -name coverage.sh`; $command" --tag cov_stats --log-dir $log_dir --queue $queue --fn-out $fn_out
	done
fi

####################################
#### call variants # the different "percent files"


perc=50
cont=20
if $local_realign; then
	echo "perform local realignment for the bam files with normal contamination $cont and percentage $perc"
	for TU_bam in `find $bamDir -name TU.wCont[0-9]0.final.RG.*$bam_ending | grep ${perc}perc | grep wCont${cont} | grep -v locRe`
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
			echo "submit local realignment"
			sleep 20
			outBam=${TU_bam%.bam}.locRealign.bam
			command="source `find $SNV_DIR -name loc_realign.sh` ; source `find $SNV_DIR -name paths.sh`; loc_realign $TU_bam $outBam $fn_genome --create-index --check-fexist "
			submit "$command" --fn-out $fn_out --tag locRealign`basename $TU_bam` --work-dir $SNV_DIR --queue $queue
		fi

		fn_out=$log_dir/loc_realign_`basename ${NO_bam%.bam}`.o
		if [ -f $fn_out ]; then
			echo "$fn_out already exists."
		else
			echo "submit local realignment"
			sleep 20
			outBamNO=${NO_bam%.bam}.locRealign.bam
			command="source `find $SNV_DIR -name loc_realign.sh`; source `find $SNV_DIR -name paths.sh`; loc_realign $NO_bam $outBamNO $fn_genome --create-index --check-fexist "
			submit "$command" --fn-out $fn_out --tag locRealign`basename $NO_bam` --work-dir $SNV_DIR --queue $queue
		fi
	done
	wait_for_jobs locRealign
fi


echo start running variant callers... 
for TU_bam in `ls $bamDir/TU.wCont[0-9]0.final.RG.*$bam_ending | grep -v locRe | grep perc`;  
do
	NO_bam=`echo $TU_bam | sed 's/TU/NO/' | sed 's/.wCont[0-9]0./_/'`;
	if [ ! -f $NO_bam ]; then
		NO_bam=$(echo $NO_bam | sed 's/final.RG./final.RG.true./') # the normal true bams have an extra "true" in their name
	fi
	echo "TU_bam = $TU_bam"
	echo "NO_bam = $NO_bam"
	perc=$( echo $TU_bam | tr '.' '\n' | grep perc | sed 's/perc//' )
	cont=$( echo $TU_bam | tr '.' '\n' | grep Cont | sed 's/wCont//')
	echo perc=$perc
	echo cont=$cont

	# for debugging purposes
	if false; then
		if [ ! $perc == 12 -o ! $cont == 20 ]; then
			continue
		fi
	fi

	if [[ $bamDir == *sn_k1* || $bamDir == *Aneuploidy* ]]; then 
		if [ ! $perc == 50 -o ! $cont == 20 ]; then
			continue
		fi
	fi

	if $local_realign; then
		if [ ! $perc == 50 -o ! $cont == 20 ]; then
			continue
		fi
		TU_bam=${TU_bam%.bam}.locRealign.bam
		NO_bam=${NO_bam%.bam}.locRealign.bam
	fi

	sample_name=TU.${perc}perc
	log_dir=$variantDir/logs/
	mkdir -p $log_dir

	#gatk
	################################################################################
	fn_out=$log_dir/gatk_`basename ${TU_bam%.bam}`_paired.o
	if [ -f $fn_out ]; then
		echo "$fn_out already exists."
		# check if error happened
		error_state=$(grep -l ERROR ${fn_out%.o}.e)
		if [ "$error_state" == "${fn_out%.o}.e" ]; then
			echo "ERROR happened with GATK for bam file $TU_bam -> submit again..."
			output=$variantDir/`basename $TU_bam | sed 's/.bam//'`.gatk_SNVs.raw.vcf
			filt_output=$variantDir/`basename $TU_bam | sed 's/.bam//'`.gatk_SNVs.raw.filt.vcf
			output_idx=${output}.idx
			filt_output_idx=${filt_output}.idx
			mv $output ${output}_ERROR
			mv $filt_output ${filt_output}_ERROR
			mv $output_idx ${output_idx}_ERROR
			mv $filt_output_idx ${filt_output_idx}_ERROR
			command="source `find $SNV_DIR -name gatkUnifiedGT_variants.sh` ; source `find $SNV_DIR -name paths.sh` ; echo $TU_bam >> $fn_out; hostname >> $fn_out ; gatkUnifiedGT_variants $TU_bam '' $fn_genome $variantDir $NO_bam"
			submit "$command" --tag gatk_$perc --fn-out $fn_out --queue $queue
		fi
	else
		command="source `find $SNV_DIR -name gatkUnifiedGT_variants.sh` ; source `find $SNV_DIR -name paths.sh` ; echo $TU_bam >> $fn_out; hostname >> $fn_out ; gatkUnifiedGT_variants $TU_bam '' $fn_genome $variantDir $NO_bam" 
		submit "$command"  --tag gatk_$perc --fn-out $fn_out --queue $queue
	fi	

	#somaticSniper
	################################################################################
	fn_out=${log_dir}somSniper_`basename ${TU_bam%.bam}`.o
	if [ -f $fn_out ]; then
			echo "$fn_out already exists."
	else
		command="source `find $SNV_DIR -name somaticSniper_variants.sh` ; source `find $SNV_DIR -name paths.sh` ; echo $TU_bam >> $fn_out; hostname >> $fn_out ; somaticSniper_variants $TU_bam '' $fn_genome $variantDir $NO_bam" 
		submit "$command" --fn-out $fn_out --tag somSniper$perc --queue $queue
	fi


	#samtools_1_2
	################################################################################
	fn_out=$log_dir/samvar_1_2_`basename ${TU_bam%.bam}`.o
	if [ -f $fn_out ]; then
			echo "$fn_out already exists."
	else
		command="source `find $SNV_DIR -name samtools_1_2_variants.sh` ; source `find $SNV_DIR -name paths.sh`; samtools_1_2_variants $TU_bam xxx $fn_genome $variantDir $NO_bam --opts --no-E "
		submit "$command" --fn-out $fn_out --tag sam_1_2$perc --work-dir $SNV_DIR --queue $queue
	fi
	#echo "wait a little bit such that the next time samtools_1_2 is invoked, it does not start again with the separate TU, NO calls"
	#sleep 10
	## also for the paired calls TU+NO together:
	#fn_out=$log_dir/samvar_1_2_pairedTUNO_`basename ${TU_bam%.bam}`.o
	#if [ -f $fn_out ]; then
	#		echo "$fn_out already exists."
	#else
	#	command="source `find $SNV_DIR -name samtools_1_2_variants.sh` ; source `find $SNV_DIR -name paths.sh`; samtools_1_2_variants $TU_bam xxx $fn_genome $variantDir $NO_bam --opts --no-E --paired"
	#	submit "$command" --fn-out $fn_out --tag sam_1_2$perc --work-dir $SNV_DIR --queue $queue
	#fi


	#varscan2
	################################################################################
	fn_out=$log_dir/varscan_`basename ${TU_bam%.bam}`.o
	if [ -f $fn_out ]; then
			echo "$fn_out already exists."
	else
		command="source `find $SNV_DIR -name varscan2_variants.sh` ; source `find $SNV_DIR -name paths.sh`; varscan2_variants $TU_bam xxx $fn_genome $variantDir $NO_bam"
		submit "$command" --fn-out $fn_out --tag varsc$perc --work-dir $SNV_DIR --queue $queue
	fi
		
	#muTect
	################################################################################
	fn_out=${log_dir}muTect_`basename ${TU_bam%.bam}`.o
	if [ -f $fn_out ]; then
		echo "$fn_out already exists."
		# check if error happened
		error_state=$(grep -l ERROR ${fn_out%.o}.e)
		if [ "$error_state" == "${fn_out%.o}.e" ]; then
			echo "ERROR happened with MuTect for bam file $TU_bam -> submit again..."
			output=$variantDir/`basename $TU_bam | sed 's/.bam//'`.muTect_SNVs_Raw.txt
			outputCoverage=$variantDir/`basename $TU_bam | sed 's/.bam//'`.muTect_coverage.wig.txt
			outputVCF=$variantDir/`basename $TU_bam | sed 's/.bam//'`.muTect_SNVs_Raw.vcf
			mv $output ${output}_ERROR
			mv $outputCoverage ${outputCoverage}_ERROR
			mv $outputVCF ${outputVCF}_ERROR

			command="source `find $SNV_DIR -name muTect_variants.sh` ; source `find $SNV_DIR -name paths.sh`; muTect_variants $TU_bam xxx $fn_genome $variantDir $NO_bam"
			submit "$command" --fn-out $fn_out --tag muTect$perc --work-dir $SNV_DIR --queue $queue
		fi
	else
		command="source `find $SNV_DIR -name muTect_variants.sh` ; source `find $SNV_DIR -name paths.sh`; muTect_variants $TU_bam xxx $fn_genome $variantDir $NO_bam"
		submit "$command" --fn-out $fn_out --tag muTect$perc --work-dir $SNV_DIR --queue $queue
	fi

	#deepSNV
	################################################################################
	fn_out=${log_dir}deepSNV_job_`basename ${TU_bam%.bam}`.o
	if [ -f $fn_out ]; then
		echo "$fn_out already exists."
	else
		#command="source `find $SNV_DIR -name deepSNV_variants.sh` ; source `find $SNV_DIR -name paths.sh` ; deepSNV_variants $TU_bam $fn_genome $variantDir $NO_bam $bed_file"
		bed_file_deep=${bed_file%.bed}_woLastColumn.bed # deepSNV needs just exaclty the three columns: chr, startPos, endPos
		command="source `find $SNV_DIR -name deepSNV_automated_variants.sh` ; source `find $SNV_DIR -name paths.sh` ; deepSNV_automated_variants $bed_file_deep $TU_bam $NO_bam $variantDir --also-vcf true"
		submit "$command" --fn-out $fn_out --tag deepSNV$perc --work-dir $SNV_DIR --queue $queue
	fi

	#jointSNVMix2
	################################################################################
	fn_out=${log_dir}jointSNVMix2_job_`basename ${TU_bam%.bam}`.o
	if [ -f $fn_out ]; then
		echo "$fn_out already exists."
	else
		command="source `find $SNV_DIR -name jointSNVMix2_variants.sh` ; source `find $SNV_DIR -name paths.sh` ; jointSNVMix2_variants $TU_bam '' $fn_genome $variantDir $NO_bam $bed_file"
		submit "$command" --fn-out $fn_out --tag jointSNVMix2$perc --work-dir $SNV_DIR --queue $queue
	fi


	#gatk HaplotypeCaller
	################################################################################
	fn_out=$log_dir/gatkHPCaller_`basename ${TU_bam%.bam}`_paired.o
 	if [ -f $fn_out ]; then
		echo "$fn_out already exists."
		## check if error happened
		error_state=$(grep -l ERROR ${fn_out%.o}.e)
		if [ "$error_state" == "${fn_out%.o}.e" ]; then
			echo "ERROR happened with GATK for bam file $TU_bam -> submit again..."
			output=$variantDir/`basename $TU_bam | sed 's/.bam//'`.gatkHPCaller_SNVs.raw.vcf
			mv $output ${output}_ERROR
			command="source `find $SNV_DIR -name gatkHPCaller_variants.sh` ; source `find $SNV_DIR -name paths.sh` ; echo $TU_bam >> $fn_out; hostname >> $fn_out ; gatkHPCaller_variants $TU_bam '' $fn_genome $variantDir $NO_bam"
			submit "$command" --tag gatkHPCaller_`basename ${TU_bam%.bam}` --fn-out $fn_out --queue $queue
		fi
	else
		command="source `find $SNV_DIR -name gatkHPCaller_variants.sh` ; source `find $SNV_DIR -name paths.sh` ; echo $TU_bam >> $fn_out; hostname >> $fn_out ; gatkHPCaller_variants $TU_bam '' $fn_genome $variantDir $NO_bam" 
		submit "$command"  --tag gatkHPCaller_`basename ${TU_bam%.bam}` --fn-out $fn_out --queue $queue
	fi

	# sinvict
	################################################################################
	fn_out=${log_dir}sinvict_job_`basename ${TU_bam%.bam}`.o
	if [ -f $fn_out ]; then
		echo "$fn_out already exists."
	else
		command="source `find $SNV_DIR -name sinvict_variants.sh`; source `find $SNV_DIR -name paths.sh` ; sinvict_variants $TU_bam $fn_genome $variantDir $NO_bam $bed_file "
		submit "$command" --fn-out $fn_out --tag sinvict$perc --work-dir $SNV_DIR --queue $queue
	fi

done
exit 0




