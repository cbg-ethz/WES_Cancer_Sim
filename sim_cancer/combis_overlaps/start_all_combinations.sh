#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`

numToolsinCombi=$1
if [ -z "$1" ]; then
	echo "Usage: $0 <numToolsinCombi> "
	exit 0
fi
echo "numToolsinCombi=$numToolsinCombi"

align_dir=$base_dir/alignments_vsn_k20/
combineScript=$currScriptDir/combininglistsbetter_Automated.R
outDir=$align_dir/variants_combined/
mkdir -p $outDir

## We take from each tool the version that performs best
joint=$align_dir/variants_default/TU.wCont20.final.RG.50perc.jointSNVMix2_SNVs_Raw.vcf
deep=$align_dir/variants_default/eval_10000_binomTest_160527/TU.wCont20.final.RG.50perc_NO_final.RG.50perc_alternative_two.sided.deepSNV.vcf_binomTestSomatic.vcf
varscan=$align_dir/variants_tuned/TU.wCont20.final.RG.50perc.bam_minvarfreq0.02_varscan2.txt.snp.Somatic_qual.vcf
gatkUG=$align_dir/variants_default/TU.wCont20.final.RG.50perc.gatk_SNVs.raw.SOMATIC.vcf
gatkHP=$align_dir/variants_default/TU.wCont20.final.RG.50perc.gatkHPCaller_SNVs.raw.SOMATIC.rewritten.vcf
somSniper=$align_dir/variants_default/TU.wCont20.final.RG.50perc.somaticSniper_SNVs_Raw_qual.noComma.vcf
SAMtoolsgzip=$align_dir/variants_tuned/TU.wCont20.final.RG.50perc.option_C200_noE_samvar_1_2.SOMATIC.vcf.gz
SAMtools=${SAMtoolsgzip%.gz}
muTect=$align_dir/variants_default/eval_10000_binomTest_160527/TU.wCont20.final.RG.50perc.muTect_SNVs_Raw.vcf_binomTestSomatic.vcf
##

# unzip samtools vcf if necessary
if [ ! -f $SAMtools ]; then
	command="gunzip $SAMtoolsgzip"
	echo $command
	eval $command
fi

declare -a all_tools=("$deep" "$varscan" "$gatkUG" "$gatkHP" "$somSniper" "$SAMtools" "$muTect" "$joint")
declare -a tools_names=('deepSNV' 'varscan2' 'gatk_SNVs' 'gatkHPCaller' 'somaticSniper' 'samvar_1_2' 'muTect' 'jointSNVMix2')
for i in `seq 1 8`; do
	echo ${all_tools[i-1]}
done

if [ $numToolsinCombi -eq 2 ]; then
	# all two-combinations
	for idxTool1 in `seq 0 7`; do
		for idxTool2 in `seq 1 7`; do
			doThisCombi=true
			tool1=${all_tools[idxTool1]}
			tool2=${all_tools[idxTool2]}

			# check that this combination does not exist yet, and that we have two different tools
			name1=${tools_names[idxTool1]}
			name2=${tools_names[idxTool2]}
			if [ "$tool1" == "$tool2" ]; then
				continue
			fi
			declare -a combinations=("$outDir/combinedlistbetter_2_${name1}_${name2}_prod.vcf" "$outDir/combinedlistbetter_2_${name2}_${name1}_prod.vcf")
			for comb_cnt in `seq 1 ${#combinations[@]}`; do
				if [ -f $combinations ]; then
					doThisCombi=false
				fi
			done
			if $doThisCombi; then
				command="$combineScript $outDir 2 $tool1 $tool2"
				echo $command
				eval $command
			fi
		done
	done
fi	

if [ $numToolsinCombi -eq 3 ]; then
	# all three-combinations
	for idxTool1 in `seq 0 7`; do
		for idxTool2 in `seq 1 7`; do
			for idxTool3 in `seq 2 7`; do
				tool1=${all_tools[idxTool1]}
				tool2=${all_tools[idxTool2]}
				tool3=${all_tools[idxTool3]}
	
				# check that this combination does not exist yet, and that we have three different tools
				declare -a curr_names=("${tools_names[idxTool1]}" "${tools_names[idxTool2]}" "${tools_names[idxTool3]}")
				if [ "$tool1" == "$tool2" -o "$tool1" == "$tool3" -o "$tool2" == "$tool3" ]; then
					continue
				fi
				for i in `seq 0 2`; do
					for j in `seq 0 2`; do
						for k in `seq 0 2`; do
							combination=$outDir/combinedlistbetter_3_${curr_names[i]}_${curr_names[j]}_${curr_names[k]}_prod.vcf
							if [ -f $combination ]; then
								continue 4 ## this continues the loop with idxTool3
							fi
						done
					done
				done
				
				# run it 
				command="$combineScript $outDir 3 $tool1 $tool2 $tool3"
				echo $command
				eval $command
			done
		done
	done
fi

if [ $numToolsinCombi -eq 4 ]; then
	# all four-combinations
	for idxTool1 in `seq 0 7`; do
		for idxTool2 in `seq 1 7`; do
			for idxTool3 in `seq 2 7`; do
				for idxTool4 in `seq 3 7`; do
					tool1=${all_tools[idxTool1]}
					tool2=${all_tools[idxTool2]}
					tool3=${all_tools[idxTool3]}
					tool4=${all_tools[idxTool4]}
	
					# check that this combination does not exist yet, and that we have only different tools
					continueLoop=false
					declare -a curr_names=("${tools_names[idxTool1]}" "${tools_names[idxTool2]}" "${tools_names[idxTool3]}" "${tools_names[idxTool4]}")
					echo "Current tools: ${curr_names[*]}"
					numTools=`echo ${curr_names[*]} | tr ' ' '\n' | sort | uniq | wc -l`
					if [ $numTools != $numToolsinCombi ]; then
						continueLoop=true
					fi
					for i in `seq 1 ${#curr_names[@]}`; do
						for j in `seq 1 ${#curr_names[@]}`; do
							for k in `seq 1 ${#curr_names[@]}`; do
								for m in `seq 1 ${#curr_names[@]}`; do
									combination=$outDir/combinedlistbetter_4_${curr_names[i-1]}_${curr_names[j-1]}_${curr_names[k-1]}_${curr_names[m-1]}_prod.vcf
									if [ -f $combination ]; then
										continueLoop=true
									fi
								done
							done
						done
					done
			
					if $continueLoop; then
						echo "continue"
						continue
					fi
	
					# run it 
					command="$combineScript $outDir 4 $tool1 $tool2 $tool3 $tool4"
					echo $command
					eval $command
				done
			done
		done
	done
fi

if [ $numToolsinCombi -eq 5 ]; then
	# all five-combinations
	for idxTool1 in `seq 0 7`; do
		for idxTool2 in `seq 1 7`; do
			for idxTool3 in `seq 2 7`; do
				for idxTool4 in `seq 3 7`; do
					for idxTool5 in `seq 4 7`; do
						tool1=${all_tools[idxTool1]}
						tool2=${all_tools[idxTool2]}
						tool3=${all_tools[idxTool3]}
						tool4=${all_tools[idxTool4]}
						tool5=${all_tools[idxTool5]}

						# check that this combination does not exist yet, and that we have only different tools
						continueLoop=false
						declare -a curr_names=("${tools_names[idxTool1]}" "${tools_names[idxTool2]}" "${tools_names[idxTool3]}" "${tools_names[idxTool4]}" "${tools_names[idxTool5]}")
						echo "Current tools: ${curr_names[*]}"

						numTools=`echo ${curr_names[*]} | tr ' ' '\n' | sort | uniq | wc -l`
						if [ $numTools != $numToolsinCombi ]; then
							continueLoop=true
						fi
						for i in `seq 1 ${#curr_names[@]}`; do
							for j in `seq 1 ${#curr_names[@]}`; do
								for k in `seq 1 ${#curr_names[@]}`; do
									for m in `seq 1 ${#curr_names[@]}`; do
										for n in `seq 1 ${#curr_names[@]}`; do
											combination=$outDir/combinedlistbetter_5_${curr_names[i-1]}_${curr_names[j-1]}_${curr_names[k-1]}_${curr_names[m-1]}_${curr_names[n-1]}_prod.vcf
	
											if [ -f $combination ]; then
												#echo "continue because combi already exists: $combination"
												continueLoop=true
											fi
										done
									done
								done
							done
						done
						if $continueLoop; then
							echo "continue"
							continue
						fi
						
						# run it 
						command="$combineScript $outDir 5 $tool1 $tool2 $tool3 $tool4 $tool5"
						echo $command
						eval $command
					done
				done
			done
		done
	done
fi


if [ $numToolsinCombi -eq 6 ]; then
	# all six-combinations - there are 28
	for numIdx in 01 02 03 04 05 06 07 12 13 14 15 16 17 23 24 25 26 27 34 35 36 37 45 46 47 56 57 67; do
	# ('deepSNV' 'varscan2' 'gatk_SNVs' 'gatkHPCaller' 'somaticSniper' 'samvar_1_2' 'muTect' 'jointSNVMix2')
	#     0         1        2               3              4             5           6          7

		idxNot1=`echo $numIdx | grep -o . | awk 'NR==1'`
		idxNot2=`echo $numIdx | grep -o . | awk 'NR==2'`
		
		cnt=0
		allToolsToCombine=""
		for allNums in 0 1 2 3 4 5 6 7; do
			if [ "$allNums" != "$idxNot1" -a "$allNums" != "$idxNot2" ]; then
				allToolsToCombine="$allToolsToCombine ${all_tools[$allNums]}"
				cnt=$((cnt+1))
			fi
		done
		# run it 
		command="$combineScript $outDir 6 ${allToolsToCombine}"
		echo $command; 	eval $command
	done
fi


if [ $numToolsinCombi -eq 7 ]; then
	# all seven-combinations - there are only 8

	# without tool 7
	tool1=${all_tools[0]}
	tool2=${all_tools[1]}
	tool3=${all_tools[2]}
	tool4=${all_tools[3]}
	tool5=${all_tools[4]}
	tool6=${all_tools[5]}
	tool7=${all_tools[6]}	
	# run it 
	command="$combineScript $outDir 7 $tool1 $tool2 $tool3 $tool4 $tool5 $tool6 $tool7"
	echo $command
	eval $command

	# without tool 6
	tool1=${all_tools[0]}
	tool2=${all_tools[1]}
	tool3=${all_tools[2]}
	tool4=${all_tools[3]}
	tool5=${all_tools[4]}
	tool6=${all_tools[5]}
	tool7=${all_tools[7]}
	# run it 
	command="$combineScript $outDir 7 $tool1 $tool2 $tool3 $tool4 $tool5 $tool6 $tool7"
	echo $command
	eval $command

	# without tool 5
	tool1=${all_tools[0]}
	tool2=${all_tools[1]}
	tool3=${all_tools[2]}
	tool4=${all_tools[3]}
	tool5=${all_tools[4]}
	tool6=${all_tools[6]}
	tool7=${all_tools[7]}
	# run it 
	command="$combineScript $outDir 7 $tool1 $tool2 $tool3 $tool4 $tool5 $tool6 $tool7"
	echo $command
	eval $command

	# without tool 4
	tool1=${all_tools[0]}
	tool2=${all_tools[1]}
	tool3=${all_tools[2]}
	tool4=${all_tools[3]}
	tool5=${all_tools[5]}
	tool6=${all_tools[6]}
	tool7=${all_tools[7]}
	# run it 
	command="$combineScript $outDir 7 $tool1 $tool2 $tool3 $tool4 $tool5 $tool6 $tool7"
	echo $command
	eval $command

	# without tool 3
	tool1=${all_tools[0]}
	tool2=${all_tools[1]}
	tool3=${all_tools[2]}
	tool4=${all_tools[4]}
	tool5=${all_tools[5]}
	tool6=${all_tools[6]}
	tool7=${all_tools[7]}
	# run it 
	command="$combineScript $outDir 7 $tool1 $tool2 $tool3 $tool4 $tool5 $tool6 $tool7"
	echo $command
	eval $command

	# without tool 2
	tool1=${all_tools[0]}
	tool2=${all_tools[1]}
	tool3=${all_tools[3]}
	tool4=${all_tools[4]}
	tool5=${all_tools[5]}
	tool6=${all_tools[6]}
	tool7=${all_tools[7]}
	# run it 
	command="$combineScript $outDir 7 $tool1 $tool2 $tool3 $tool4 $tool5 $tool6 $tool7"
	echo $command
	eval $command

	# without tool 1
	tool1=${all_tools[0]}
	tool2=${all_tools[2]}
	tool3=${all_tools[3]}
	tool4=${all_tools[4]}
	tool5=${all_tools[5]}
	tool6=${all_tools[6]}
	tool7=${all_tools[7]}
	# run it 
	command="$combineScript $outDir 7 $tool1 $tool2 $tool3 $tool4 $tool5 $tool6 $tool7"
	echo $command
	eval $command

	# without tool 0
	tool1=${all_tools[1]}
	tool2=${all_tools[2]}
	tool3=${all_tools[3]}
	tool4=${all_tools[4]}
	tool5=${all_tools[5]}
	tool6=${all_tools[6]}
	tool7=${all_tools[7]}
	# run it 
	command="$combineScript $outDir 7 $tool1 $tool2 $tool3 $tool4 $tool5 $tool6 $tool7"
	echo $command
	eval $command
fi


if [ $numToolsinCombi -eq 8 ]; then
	# all eight tools
	command="$combineScript $outDir 8 ${all_tools[*]}"
	echo $command
	eval $command
fi


