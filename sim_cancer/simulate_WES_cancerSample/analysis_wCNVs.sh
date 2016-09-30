#!/bin/bash

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`
source `find ${gitDir} -name paths.sh`

# generate the overlaps between all subclones to obtain the intermediate clones vcf:
./generate_vcf_files_allSubclones.sh $vcf_sim_dir $pipeline_eval/final_vcfs_0.25_subclones/ 
# (uses the program: SNV_Overlapper_Ordered_PlusAlleles.java)

# run control-freec to obtain cnv calls from the true cancer file; $base_dir_wCNVs/../generate_vcfs/ contains the bam files 
./run_control_freec_for_cnvs.sh $pipeline_eval/generate_vcfs/  $pipeline_eval/CNVs/controlFreec/

# distribute the CNVs onto the tree: randomly deciding tree level, tree node, and zygosity; all children inherit the CNVs from their parents
outDir=$pipeline_eval/final_vcfs_0.25_subclones/CNVs
mkdir -p $outDir
cnvFile=$pipeline_eval/CNVs/controlFreec/19_4845-M_TU1_C3_mm_smmr.rmdup.bam_CNVs
python distribute_cnvs_on_tree.py $cnvFile $outDir

# delete all cnv lists which are empty:
cnvDir=$pipeline_eval/final_vcfs_0.25_subclones/CNVs/
for bedFile in $cnvDir/*bed; do
	if  [  ! -s $bedFile ]; then
		rm $bedFile
	fi
done

# extract all SNVs which are not in a loss for each node: (CN losses)
./extract_variants_in_cnv_regions_losses.sh $pipeline_eval/final_vcfs_0.25_subclones/CNVs/ $pipeline_eval/final_vcfs_0.25_subclones/

# additionally, loss of a whole chromosome in a node:
# simulate loss of chromosome 6 in node 10_0:
cnvListDir=$pipeline_eval/final_vcfs_0.25_subclones/CNVs/ 
intersectedVarDir=$cnvListDir/intersectedVariants/
outDir=$intersectedVarDir/chrLoss/
mkdir -p $outDir
nodeInTree=10_0_TU
lossFile10_0_TU=$cnvListDir/cnv_losses_10_0_TU_woCHR6.bed
echo -e "chr6\t1\t171115067\tloss" > $lossFile10_0_TU ## whole chromosome 6
$intersectBed -a $intersectedVarDir/node_10_0_TU_withoutLossRegions.txt -b $lossFile10_0_TU -v > $outDir/node_${nodeInTree}_withoutLossRegions.txt
for varFile in $intersectedVarDir/*.txt; do 
	if [ ! -f $outDir/`basename $varFile` ]; then
		echo "ln -s $varFile $outDir/`basename $varFile`"
		ln -s $varFile $outDir/`basename $varFile`
	fi
done

# gain:
./extract_variants_in_cnv_regions_gains.sh $pipeline_eval/final_vcfs_0.25_subclones/CNVs/ $pipeline_eval/final_vcfs_0.25_subclones/CNVs/intersectedVariants/chrLoss/

# create final mutation lists for each subclone - with all losses and gains:
inclGainDir=$cnvListDir/intersectedVariants/chrLoss/inclGains/
outDir=$inclGainDir/listPerSubclone/
mkdir -p $outDir
for list in $inclGainDir/*[nU].txt; do
	ln -s $list $outDir/`basename ${list}`
done
for list in $inclGainDir/../*.txt; do
	subcloneName=`basename $list | sed 's/_withoutLossRegions.txt//' | sed 's/.txt//'`
	if [ ! -f $outDir/${subcloneName}.txt ]; then
		ln -s $list $outDir/${subcloneName}.txt
	fi
done

# add up all mutations for each leaf subclone, get all bed files, also for gains
cnvDir=$pipeline_eval/final_vcfs_0.25_subclones/CNVs/ 
./combine_and_sort_subclones_variants_n_regions.sh $cnvDir $cnvDir/intersectedVariants/chrLoss/inclGains/listPerSubclone/

# create fasta from vcf
vcfDir=$vcf_sim_dir_wCNVs # = $pipeline_eval/final_vcfs_0.25_subclones/CNVs/intersectedVariants/chrLoss/inclGains/listPerSubclone/finalListsAndRegions/
outDir=$vcfDir
mkdir -p $outDir
for fn_vcf in $vcfDir/*.vcf
do
	fn_ref_out=${fn_vcf%.vcf}.fa
	fn_map=${fn_vcf%.vcf}_mapping.csv
	command="./create_ref_from_vcf $fn_genome $fn_vcf $fn_ref_out $fn_map > $outDir/`basename ${fn_vcf%.vcf}`.log 2>&1 "
	echo $command
	eval $command
done

# link the normal fastqs and mapping.csv files from the original simulation: They do not have CNVs anyways:
for copy in 0 1; do
	ln -s $vcf_sim_dir/${copy}_NO_final.fa $vcf_sim_dir_wCNVs/${copy}_NO_final.fa
	ln -s $vcf_sim_dir/${copy}_NO_final_mapping.csv $vcf_sim_dir_wCNVs/${copy}_NO_final_mapping.csv
done

# generate reads from the fastas
./run_wessim_wCNVs.sh $vcf_sim_dir_wCNVs  $base_dir_wCNVs

