###
# order of scripts to run for simulation:
prepare_input_bams.sh: creates the two bam files from each of the tumor and normal bam files, one with all multimappers, and one without multimappers

create_cancer_vcfs_as_input.sh:	runs freebayes on all bam files

intersect_extract_merge_vcfs.sh: from the final tumor and normal vcf files, generate the intersections as desired; this script calls compute_overlap.cpp

run_create_cancer_vcfs.sh: inserts the variants into the tree and then, when you have a vcf for each clone, insert them into the reference genome, for each clone; this script calls create_cancer_vcfs.py and create_ref_from_vcf.cpp

run_wessim_all.sh: for each of the cancer clone genomes, it generates reads and aligns them; this script calls run_wessim.sh, process_alignments.sh, process_true_alignments.sh; first: make sure that 'select_true_alignment.cpp' in the folder ../eval_vcf/ is compiled


###
# Now, the variant callers can be started. see ../var_callers/

###
# these are additional scripts which can be run:

check_coverage.sh: for coverage statistics and qualimap

repeated_subsampling.sh: to assert the variability of the performance, repeated subsampling of the whole (100%) bam file


## revisions:

## generate the overlaps between all subclones to obtain the intermediate clones vcf:
## ./generate_vcf_files_allSubclones.sh  /home/arhofman/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25 ~/scratch/arhofman/exome_seq_pipeline_eval/final_vcfs_0.25_subclones/
## uses the program: SNV_Overlapper_Ordered_PlusAlleles.java

## run control-freec to obtain cnv calls from the true cancer file
## ./run_control_freec_for_cnvs.sh /links/grid/scratch/arhofman/exome_seq_pipeline_eval/generate_vcfs/ /links/grid/scratch/arhofman/exome_seq_pipeline_eval/CNVs/controlFreec/

## make sure that there are now overlapping CNVs (should be the case anyway)
## awk 'BEGIN{chr="";end=0}{if(chr==$1 && $2<end){print};chr=$1;end=$3}' 19_4845-M_TU1_C3_mm_smmr.rmdup.bam_CNVs

## distribute the CNVs onto the tree: randomly deciding tree level, tree node, and zygosity; all children inherit the CNVs from their parents
outDir=~/scratch/arhofman/exome_seq_pipeline_eval/final_vcfs_0.25_subclones/CNVs
mkdir -p $outDir
cnvFile=/links/grid/scratch/arhofman/exome_seq_pipeline_eval/CNVs/controlFreec/19_4845-M_TU1_C3_mm_smmr.rmdup.bam_CNVs
python distribute_cnvs_on_tree.py $cnvFile $outDir

## delete all cnv lists which are empty:
cnvDir=~/scratch/arhofman/exome_seq_pipeline_eval/final_vcfs_0.25_subclones/CNVs/
for bedFile in $cnvDir/*bed; do
	if  [  ! -s $bedFile ]; then
		rm $bedFile
	fi
done

## extract all SNVs which are not in a loss for each node: (CN losses)
./extract_variants_in_cnv_regions_losses.sh /home/arhofman/scratch/arhofman/exome_seq_pipeline_eval/final_vcfs_0.25_subclones/CNVs/  /home/arhofman/scratch/arhofman/exome_seq_pipeline_eval/final_vcfs_0.25_subclones/

## additionally, loss of a whole chromosome in a node:
# simulate loss of chromosome 6 in node 10_0:
cnvListDir=/home/arhofman/scratch/arhofman/exome_seq_pipeline_eval/final_vcfs_0.25_subclones/CNVs/
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

## gain:
./extract_variants_in_cnv_regions_gains.sh /home/arhofman/scratch/arhofman/exome_seq_pipeline_eval/final_vcfs_0.25_subclones/CNVs/ /home/arhofman/scratch/arhofman/exome_seq_pipeline_eval/final_vcfs_0.25_subclones/CNVs/intersectedVariants/chrLoss/

## create final mutation lists for each subclone - with all losses and gains:
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


## add up all mutations for each leaf subclone, get all bed files, also for gains
cnvDir=~/scratch/arhofman/exome_seq_pipeline_eval/final_vcfs_0.25_subclones/CNVs/
./combine_and_sort_subclones_variants_n_regions.sh $cnvDir $cnvDir/intersectedVariants/chrLoss/inclGains/listPerSubclone/


## create fasta from vcf
source ../../paths.sh
vcfDir=/home/arhofman/scratch/arhofman/exome_seq_pipeline_eval/final_vcfs_0.25_subclones/CNVs/intersectedVariants/chrLoss/inclGains/listPerSubclone/finalListsAndRegions/
outDir=$vcfDir
mkdir -p $outDir
create_ref_from_vcfFolder=~/shared/singerj/git/sim_cancer_paper/
cd $create_ref_from_vcfFolder
	for fn_vcf in $vcfDir/*.vcf
	do
		fn_ref_out=${fn_vcf%.vcf}.fa
		fn_map=${fn_vcf%.vcf}_mapping.csv
		command="./create_ref_from_vcf $fn_genome $fn_vcf $fn_ref_out $fn_map > $outDir/`basename ${fn_vcf%.vcf}`.log 2>&1 "
		echo $command
		eval $command
	done
cd -

## link the normal fastqs and mapping.csv files from the original simulation: They do not have CNVs anyways:
for copy in 0 1; do
	ln -s /home/arhofman/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/${copy}_NO_final.fa /net/bs-gridfs/gridexport/scratch/arhofman/exome_seq_pipeline_eval/final_vcfs_0.25_subclones/CNVs/intersectedVariants/chrLoss/inclGains/listPerSubclone/finalListsAndRegions/${copy}_NO_final.fa
	ln -s /home/arhofman/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/${copy}_NO_final_mapping.csv /net/bs-gridfs/gridexport/scratch/arhofman/exome_seq_pipeline_eval/final_vcfs_0.25_subclones/CNVs/intersectedVariants/chrLoss/inclGains/listPerSubclone/finalListsAndRegions/${copy}_NO_final_mapping.csv
done

## generate reads from the fastas
./run_wessim_wCNVs.sh /net/bs-gridfs/gridexport/scratch/arhofman/exome_seq_pipeline_eval/final_vcfs_0.25_subclones/CNVs/intersectedVariants/chrLoss/inclGains/listPerSubclone/finalListsAndRegions/  /net/bs-gridfs/gridexport/scratch/arhofman/exome_seq_pipeline_eval/simulation_wCNVs_n_Aneuploidy/



