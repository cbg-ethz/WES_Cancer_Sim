#!/bin/bash


currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`


simDir=$base_dir/alignments_vsn_k20/
simQualimap=$simDir/qualimap
declare -a coverageProfile=( "$simQualimap/TU.wCont20.final.RG.50perc/TU.wCont20.final.RG.50perc.txt")

for fileIDX in `seq 1 ${#coverageProfile[@]}`; do
	myCovProfile=${coverageProfile[fileIDX-1]}
	echo $myCovProfile
	meanCov=$(cat $myCovProfile | awk 'NR==71' | awk '{print $4}')
	echo $meanCov
	Rsummary=${myCovProfile%.txt}_inputR.txt
	if [ ! -f $Rsummary ]; then
		cat $myCovProfile | awk 'NR>=74' | awk 'NR <=51' | awk '{print $11"\t"$4}' | sed 's/%//g' | sed 's/X//g' > $Rsummary
	fi
	sampleName=`basename $myCovProfile`
	$RscriptPath $gitDir/sim_cancer/eval_vcf/plot_coverage_profiles.R $Rsummary $meanCov $sampleName
done


covProfileTotal=$simDir/quality_control_default/TU.wCont20.final.RG.50perc.genome_frac_cov.txt
meanCov=$(cat $covProfileTotal | awk 'NR==3' | awk '{print $1}')
medianCov=$(cat $covProfileTotal | awk 'NR==4' | awk '{print $1}')
Rsummary=${covProfileTotal%.txt}_inputR.txt
if [ ! -f $Rsummary ]; then
	cat $covProfileTotal | awk 'NR>=7' | awk '{print $1"\t"$2}' > $Rsummary
fi
sampleName=`basename $covProfileTotal`
$RscriptPath $gitDir/sim_cancer/eval_vcf/plot_coverage_profiles.R $Rsummary $meanCov $sampleName $medianCov


