#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../ -name paths.sh`

numToolsinCombi=$1
if [ -z "$1" ]; then
	echo "Usage: $0 <numToolsinCombi> "
	exit 0
fi
echo "numToolsinCombi=$numToolsinCombi"

align_dir=$base_dir/alignments_vsn_k20/
combineScript=$currScriptDir/combininglistsbetter_Automated.R
#outDir=$align_dir/variants_combined/ 
outDir=$align_dir/variants_combined_sinvict/
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
sinvict=${align_dir}/variants_tuned/sinvictqscorecutoff60_TU.wCont20.final.RG.50perc_vs_NO_final.RG.50perc_somatic.vcf
##

# unzip samtools vcf if necessary
if [ ! -f $SAMtools ]; then
	command="gunzip $SAMtoolsgzip"
	echo $command
	eval $command
fi

declare -a all_tools=("$deep" "$varscan" "$gatkUG" "$gatkHP" "$somSniper" "$SAMtools" "$muTect" "$joint" "$sinvict")
declare -a tools_names=('deepSNV' 'varscan2' 'gatk_SNVs' 'gatkHPCaller' 'somaticSniper' 'samvar_1_2' 'muTect' 'jointSNVMix2' 'sinvict')
for i in `seq 1 9`; do
	echo ${all_tools[i-1]}
done


if [ $numToolsinCombi -eq 2 ]; then
	# all two-combinations - there are 36
	# for i in `seq 0 8`; do oneAfter=$((i+1)); if [ $oneAfter -le 8 ]; then for j in `seq $oneAfter 8`; do if [ $i == $j ]; then echo error; fi; echo "$i$j"; done;fi;done | tr '\n' ' '
	for numIdx in 01 02 03 04 05 06 07 08 12 13 14 15 16 17 18 23 24 25 26 27 28 34 35 36 37 38 45 46 47 48 56 57 58 67 68 78; do
		# ('deepSNV' 'varscan2' 'gatk_SNVs' 'gatkHPCaller' 'somaticSniper' 'samvar_1_2' 'muTect' 'jointSNVMix2' 'sinvict')
                #     0         1        2               3              4             5           6          7           8
		idxIn1=`echo $numIdx | grep -o . | awk 'NR==1'`
                idxIn2=`echo $numIdx | grep -o . | awk 'NR==2'`
		
		cnt=0
                allToolsToCombine=""
		for allNums in 0 1 2 3 4 5 6 7 8; do
                        if [ "$allNums" == "$idxIn1" -o "$allNums" == "$idxIn2" ]; then
				allToolsToCombine="$allToolsToCombine ${all_tools[$allNums]}"
				cnt=$((cnt+1))
			fi
		done
		# run it
                command="$combineScript $outDir 2 ${allToolsToCombine}"
		echo $command;  eval $command
	done
fi


if [ $numToolsinCombi -eq 3 ]; then
	# all three-combinations - there are 84
	# for i in `seq 0 8`; do oneAfter=$((i+1)); if [ $oneAfter -le 8 ]; then for j in `seq $oneAfter 8`; do if [ $i == $j ]; then echo error; fi; nextOne=$((j+1)); if [ $nextOne -le 8 ]; then for k in `seq $nextOne 8`; do if [ $k == $j ]; then echo error; fi;  echo "$i$j$k"; done; fi ; done; fi;done  | tr '\n' ' '
	for numIdx in 012 013 014 015 016 017 018 023 024 025 026 027 028 034 035 036 037 038 045 046 047 048 056 057 058 067 068 078 123 124 125 126 127 128 134 135 136 137 138 145 146 147 148 156 157 158 167 168 178 234 235 236 237 238 245 246 247 248 256 257 258 267 268 278 345 346 347 348 356 357 358 367 368 378 456 457 458 467 468 478 567 568 578 678; do
		# ('deepSNV' 'varscan2' 'gatk_SNVs' 'gatkHPCaller' 'somaticSniper' 'samvar_1_2' 'muTect' 'jointSNVMix2' 'sinvict')
		#     0         1        2               3              4             5           6          7		 8

		idxIn1=`echo $numIdx | grep -o . | awk 'NR==1'`
		idxIn2=`echo $numIdx | grep -o . | awk 'NR==2'`
		idxIn3=`echo $numIdx | grep -o . | awk 'NR==3'`
		
		cnt=0
		allToolsToCombine=""
		for allNums in 0 1 2 3 4 5 6 7 8; do
			if [ "$allNums" == "$idxIn1" -o "$allNums" == "$idxIn2" -o "$allNums" == "$idxIn3" ]; then
				allToolsToCombine="$allToolsToCombine ${all_tools[$allNums]}"
				cnt=$((cnt+1))
			fi
		done
		# run it 
		command="$combineScript $outDir 3 ${allToolsToCombine}"
		echo $command; 	eval $command
	done
fi


if [ $numToolsinCombi -eq 4 ]; then
	# all four-combinations - there are 126
	#for i in `seq 0 8`; do oneAfter=$((i+1)); if [ $oneAfter -le 8 ]; then for j in `seq $oneAfter 8`; do if [ $i == $j ]; then echo error; fi; nextOne=$((j+1)); if [ $nextOne -le 8 ]; then for k in `seq $nextOne 8`; do if [ $k == $j ]; then echo error; fi;  fourthOne=$((k+1)); if [ $fourthOne -le 8 ]; then for m in `seq $fourthOne 8`; do if [ $m == $k ]; then echo error; fi ; echo "$i$j$k$m";done; fi ; done; fi ; done; fi;done | tr '\n' ' '
	for numIdx in 0123 0124 0125 0126 0127 0128 0134 0135 0136 0137 0138 0145 0146 0147 0148 0156 0157 0158 0167 0168 0178 0234 0235 0236 0237 0238 0245 0246 0247 0248 0256 0257 0258 0267 0268 0278 0345 0346 0347 0348 0356 0357 0358 0367 0368 0378 0456 0457 0458 0467 0468 0478 0567 0568 0578 0678 1234 1235 1236 1237 1238 1245 1246 1247 1248 1256 1257 1258 1267 1268 1278 1345 1346 1347 1348 1356 1357 1358 1367 1368 1378 1456 1457 1458 1467 1468 1478 1567 1568 1578 1678 2345 2346 2347 2348 2356 2357 2358 2367 2368 2378 2456 2457 2458 2467 2468 2478 2567 2568 2578 2678 3456 3457 3458 3467 3468 3478 3567 3568 3578 3678 4567 4568 4578 4678 5678; do
		# ('deepSNV' 'varscan2' 'gatk_SNVs' 'gatkHPCaller' 'somaticSniper' 'samvar_1_2' 'muTect' 'jointSNVMix2' 'sinvict')
		#     0         1        2               3              4             5           6          7           8
	
		idxIn1=`echo $numIdx | grep -o . | awk 'NR==1'`
		idxIn2=`echo $numIdx | grep -o . | awk 'NR==2'`
		idxIn3=`echo $numIdx | grep -o . | awk 'NR==3'`
		idxIn4=`echo $numIdx | grep -o . | awk 'NR==4'`

		cnt=0
		allToolsToCombine=""
		for allNums in 0 1 2 3 4 5 6 7 8; do
			if [ "$allNums" == "$idxIn1" -o "$allNums" == "$idxIn2" -o "$allNums" == "$idxIn3" -o "$allNums" == "$idxIn4" ]; then
				allToolsToCombine="$allToolsToCombine ${all_tools[$allNums]}"
				cnt=$((cnt+1))
			fi
		done
		# run it
		command="$combineScript $outDir 4 ${allToolsToCombine}"
		echo $command;  eval $command
	done
fi



if [ $numToolsinCombi -eq 5 ]; then
	# all five-combinations - there are 126
	#for i in `seq 0 8`; do oneAfter=$((i+1)); if [ $oneAfter -le 8 ]; then for j in `seq $oneAfter 8`; do if [ $i == $j ]; then echo error; fi; nextOne=$((j+1)); if [ $nextOne -le 8 ]; then for k in `seq $nextOne 8`; do if [ $k == $j ]; then echo error; fi;  fourthOne=$((k+1)); if [ $fourthOne -le 8 ]; then for m in `seq $fourthOne 8`; do if [ $m == $k ]; then echo error; fi ; echo "$i$j$k$m";done; fi ; done; fi ; done; fi;done | tr '\n' ' '
	for numIdx in 0123 0124 0125 0126 0127 0128 0134 0135 0136 0137 0138 0145 0146 0147 0148 0156 0157 0158 0167 0168 0178 0234 0235 0236 0237 0238 0245 0246 0247 0248 0256 0257 0258 0267 0268 0278 0345 0346 0347 0348 0356 0357 0358 0367 0368 0378 0456 0457 0458 0467 0468 0478 0567 0568 0578 0678 1234 1235 1236 1237 1238 1245 1246 1247 1248 1256 1257 1258 1267 1268 1278 1345 1346 1347 1348 1356 1357 1358 1367 1368 1378 1456 1457 1458 1467 1468 1478 1567 1568 1578 1678 2345 2346 2347 2348 2356 2357 2358 2367 2368 2378 2456 2457 2458 2467 2468 2478 2567 2568 2578 2678 3456 3457 3458 3467 3468 3478 3567 3568 3578 3678 4567 4568 4578 4678 5678; do
		# ('deepSNV' 'varscan2' 'gatk_SNVs' 'gatkHPCaller' 'somaticSniper' 'samvar_1_2' 'muTect' 'jointSNVMix2' 'sinvict')
		#     0         1        2               3              4             5           6          7           8
	
		idxNot1=`echo $numIdx | grep -o . | awk 'NR==1'`
		idxNot2=`echo $numIdx | grep -o . | awk 'NR==2'`
		idxNot3=`echo $numIdx | grep -o . | awk 'NR==3'`
		idxNot4=`echo $numIdx | grep -o . | awk 'NR==4'`

		cnt=0
		allToolsToCombine=""
		for allNums in 0 1 2 3 4 5 6 7 8; do
			if [ "$allNums" != "$idxNot1" -a "$allNums" != "$idxNot2" -a "$allNums" != "$idxNot3" -a "$allNums" != "$idxNot4" ]; then
				allToolsToCombine="$allToolsToCombine ${all_tools[$allNums]}"
				cnt=$((cnt+1))
			fi
		done
		# run it
		command="$combineScript $outDir 5 ${allToolsToCombine}"
		echo $command;  eval $command
	done
fi

if [ $numToolsinCombi -eq 6 ]; then
	# all six-combinations - there are 84
	# for i in `seq 0 8`; do oneAfter=$((i+1)); if [ $oneAfter -le 8 ]; then for j in `seq $oneAfter 8`; do if [ $i == $j ]; then echo error; fi; nextOne=$((j+1)); if [ $nextOne -le 8 ]; then for k in `seq $nextOne 8`; do if [ $k == $j ]; then echo error; fi;  echo "$i$j$k"; done; fi ; done; fi;done  | tr '\n' ' '
	for numIdx in 012 013 014 015 016 017 018 023 024 025 026 027 028 034 035 036 037 038 045 046 047 048 056 057 058 067 068 078 123 124 125 126 127 128 134 135 136 137 138 145 146 147 148 156 157 158 167 168 178 234 235 236 237 238 245 246 247 248 256 257 258 267 268 278 345 346 347 348 356 357 358 367 368 378 456 457 458 467 468 478 567 568 578 678; do
		# ('deepSNV' 'varscan2' 'gatk_SNVs' 'gatkHPCaller' 'somaticSniper' 'samvar_1_2' 'muTect' 'jointSNVMix2' 'sinvict')
		#     0         1        2               3              4             5           6          7		 8

		idxNot1=`echo $numIdx | grep -o . | awk 'NR==1'`
		idxNot2=`echo $numIdx | grep -o . | awk 'NR==2'`
		idxNot3=`echo $numIdx | grep -o . | awk 'NR==3'`
		
		cnt=0
		allToolsToCombine=""
		for allNums in 0 1 2 3 4 5 6 7 8; do
			if [ "$allNums" != "$idxNot1" -a "$allNums" != "$idxNot2" -a "$allNums" != "$idxNot3" ]; then
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
	# all seven-combinations - there are 36
	# for i in `seq 0 8`; do oneAfter=$((i+1)); if [ $oneAfter -le 8 ]; then for j in `seq $oneAfter 8`; do if [ $i == $j ]; then echo error; fi; echo "$i$j"; done;fi;done | tr '\n' ' '
	for numIdx in 01 02 03 04 05 06 07 08 12 13 14 15 16 17 18 23 24 25 26 27 28 34 35 36 37 38 45 46 47 48 56 57 58 67 68 78; do
		# ('deepSNV' 'varscan2' 'gatk_SNVs' 'gatkHPCaller' 'somaticSniper' 'samvar_1_2' 'muTect' 'jointSNVMix2' 'sinvict')
                #     0         1        2               3              4             5           6          7           8
		idxNot1=`echo $numIdx | grep -o . | awk 'NR==1'`
                idxNot2=`echo $numIdx | grep -o . | awk 'NR==2'`
		
		cnt=0
                allToolsToCombine=""
		for allNums in 0 1 2 3 4 5 6 7 8; do
                        if [ "$allNums" != "$idxNot1" -a "$allNums" != "$idxNot2" ]; then
				allToolsToCombine="$allToolsToCombine ${all_tools[$allNums]}"
				cnt=$((cnt+1))
			fi
		done
		# run it
                command="$combineScript $outDir 7 ${allToolsToCombine}"
		echo $command;  eval $command
	done
fi



if [ $numToolsinCombi -eq 8 ]; then
	# all eight combinations - there are 9
	for numIdx in 0 1 2 3 4 5 6 7 8; do
		# ('deepSNV' 'varscan2' 'gatk_SNVs' 'gatkHPCaller' 'somaticSniper' 'samvar_1_2' 'muTect' 'jointSNVMix2' 'sinvict')
                #     0         1        2               3              4             5           6          7           8
		cnt=0
                allToolsToCombine=""
		for allNums in 0 1 2 3 4 5 6 7 8; do
			if [ "$allNums" != "$numIdx" ]; then
				allToolsToCombine="$allToolsToCombine ${all_tools[$allNums]}"
				cnt=$((cnt+1))
			fi
		done
		# run it
		command="$combineScript $outDir 8 ${allToolsToCombine}"
		echo $command
		eval $command
	done
fi


if [ $numToolsinCombi -eq 9 ]; then
	# all nine tools
	command="$combineScript $outDir 9 ${all_tools[*]}"
	echo $command
	eval $command
fi


