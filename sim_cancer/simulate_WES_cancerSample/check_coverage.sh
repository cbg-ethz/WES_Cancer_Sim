#/bin/bash -l

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`

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
outputDir=$2
bam_ending=$3
if [ -z "$3" ]; then
	bam_ending=perc.bam
fi

## perform quality control with qualimap
qualityOut=$outputDir/qualimap
logDir=$qualityOut/logs/
mkdir -p $qualityOut
mkdir -p $logDir

bamShort=TU.wCont20.final.RG.50perc
genomeFractionCoverage=$outputDir/quality_control_default/${bamShort}.genome_frac_cov.txt
if [ ! -f $genomeFractionCoverage ]; then
	covStatScript=`find ../ -name coverageR_quality_stats.R`
	logDir=$outputDir/quality_control_default/logs/
	fn_out=${logDir}genomeFrac_job_${bamShort}.o
	if [ ! -f $fn_out ]; then
		RStatCov=$outputDir/quality_control_default/coverage/perBaseCoverage_${bamShort}_R.txt
		if [ -f $RStatCov ]; then
			command="$covStatScript $RStatCov > $genomeFractionCoverage"
			submit "$command" --tag genomeFrac_${bamShort} --fn-out $fn_out --queue mpi01.q
		fi
	fi
fi


for curr_bam in $bamDir/TU.wCont[0-9]0.final.RG.*$bam_ending $bamDir/NO_final.RG.*$bam_ending
do
	echo $curr_bam
	fn_out=${logDir}cov_job_`basename ${curr_bam%.bam}`.o
	if [ -f $fn_out ]; then
		continue;
	fi
	currOutDir=$qualityOut/`basename ${curr_bam%.bam}`
	if [ ! -d $currOutDir ]; then
		command="$qualimap bamqc -bam $curr_bam --feature-file $bed_file_qualimap -outdir $currOutDir -outfile `basename ${curr_bam%.bam}`.report.pdf -outformat PDF -nt 1 --java-mem-size=4G"
		echo $command
		submit "$command" --tag qualimap`basename ${curr_bam%.bam}` --log-dir $logDir --queue mpi01.q --fn-out $fn_out 
		#eval $command
		#exit 0
	fi
done
exit 0
covSummary=$qualityOut/coverage_summary.txt
if [ ! -f $covSummary -o ! -s $covSummary ]; then # if does not exist or if empty
	for i in `find $qualityOut/ -name genome_results.txt`; do
		sampleName=$(echo `dirname $i | tr '/' ' ' | awk '{print $NF}'`)
		ln -s $i `dirname $i`/${sampleName}.txt 
		meanCov=$(cat $i | grep mean | grep coverageData | awk '{print $4}' | sed 's/X//')
		echo -e "$sampleName\t$meanCov" >> $covSummary
	done
fi


