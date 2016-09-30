###
# scripts to run for variant calling:
start_callers.sh: just prints out the commands to run for variant calling in default mode, variant calling with parameter tuning, repeated subsampling, local realignment, checking the coverage, ... ; prints the run commands for run_variant_callers.sh, tune_variant_callers.sh

###
# scripts to run for filtering and re-formatting of the output
filterOutGermlineComma_somSniper.sh: filters out germline mutations from the somaticSniper output; calls the program Extract_Somatic_Mult_Alleles.java

add_qual_in_VCF_for_ordering.sh: reformats the output of VarScan2 and somaticSniper. It writes the quality score of the variant into the 6th column, like is the case for a usual vcf file; calls the programs RewriteVarScanVCF.jar and RewriteSomaticSniperQualField.java

filterOutGermline_GATK_SAMtools.sh: filters out the germline mutations from the vcfs of GATK and SAMtools; calls the program WriteDifference_TU_NO_asSomaticFilter.java 

format_HPCaller.sh: re-formats the GATK HP vcf, because sometimes the column with the "GT:AD:DP:GQ:PL" details is malformed; calls the program FormatHPCallerGATK_TUonly.java

###
# Now, the vcfs of the variant callers can be evaluated. see ../eval_vcf/

# or combined/intersected. see ../combis_overlaps/


#################
# revisions:

# test sinvict:
./test_sinvict.sh ~/scratch/arhofman/exome_seq_pipeline_eval/generate_vcfs/  ~/scratch/arhofman/exome_seq_pipeline_eval/testSinvict/

# run the variant callers such that sinvict is done:
## all eight vsn k20 bams:
5043  2016-09-12 19:11:50 ./run_variant_callers.sh /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_vsn_k20/bam /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_vsn_k20
## one vsn k20 bam with loc realign true:
 5045  2016-09-12 19:12:11 ./run_variant_callers.sh /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_vsn_k20/bam /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_vsn_k20
## all eight true vsn k20 bams:
5081  2016-09-13 17:38:43 ./run_variant_callers.sh /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_vsn_k20/bam /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_vsn_k20 perc.true.bam
## the ten repeated subsampling bams: (change tag!!)
5096  2016-09-13 17:45:48 ./run_variant_callers.sh /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_vsn_k20/bam/repeatedSubsampling/ /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_vsn_k20/bam/  perc.[0-9]\*.bam

## one sn k1 bam file (50 perc)
5112  2016-09-13 18:18:00 ./run_variant_callers.sh /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_sn_k1/bam /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_sn_k1
## one true sn k1 bam file (50 perc)
 5115  2016-09-13 18:18:40 ./run_variant_callers.sh /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_sn_k1/bam /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_sn_k1 perc.true.bam

--> 29 runs

# plus 4 runs: tuning with respect to qscore cutoff:
./tune_variant_callers.sh /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_vsn_k20/bam /links/grid/scratch/singerj/sim_cancer_exome/data/subsample/final_vcfs_0.25/tree_blat_i70_s60/alignments_vsn_k20

## start all 9 variant callers on the newly simulated data
## all 9 callers for the vsn k20 50 perc 20 cont file
./run_variant_callers.sh /net/bs-gridfs/gridexport/scratch/arhofman/exome_seq_pipeline_eval/simulation_wCNVs_n_Aneuploidy/alignments_vsn_k20/bam/ /net/bs-gridfs/gridexport/scratch/arhofman/exome_seq_pipeline_eval/simulation_wCNVs_n_Aneuploidy/alignments_vsn_k20

# all 9 callers for the true bam the vsn k20 50 perc 20 cont
./run_variant_callers.sh /net/bs-gridfs/gridexport/scratch/arhofman/exome_seq_pipeline_eval/simulation_wCNVs_n_Aneuploidy/alignments_vsn_k20/bam/ /net/bs-gridfs/gridexport/scratch/arhofman/exome_seq_pipeline_eval/simulation_wCNVs_n_Aneuploidy/alignments_vsn_k20 perc.true.bam


## rewrite/filter the new variant calls:
# somaticsniper and varscan: quality field
source ../../paths.sh
for somaticSniperVCF in `find $base_dir_wCNVs/alignments_*/variants_default/ -name *somaticSniper_SNVs_Raw.vcf`
do
	if [ ! -f ${somaticSniperVCF%.vcf}_qual.vcf ]; then
		echo "$somaticSniperVCF"

		command="java RewriteSomaticSniperQualField $somaticSniperVCF > ${somaticSniperVCF%.vcf}_qual.vcf "
		echo $command
		eval $command
	fi
done
for varscanVCF in `find $base_dir_wCNVs/alignments_*/variants*/ -name *varscan2.txt.snp.Somatic.vcf`
do
	if [ ! -f ${varscanVCF%.vcf}_qual.vcf ]; then
		echo "VarScan2"
		command="java -jar RewriteVarScanVCF.jar $varscanVCF > ${varscanVCF%.vcf}_qual.vcf"
		echo $command
		eval $command
	fi
done


## gatk and samtools: filter out germline:
for SAMVCF in `find $base_dir_wCNVs/alignments_vsn_k20/variants_default/ -name TU*samvar*.vcf.gz | grep -v eval_1000 | grep -v GERMLINE | grep -v SOMATIC | grep -v combined | grep -v SOMATIC_paired`; do
        NoVCF=${SAMVCF%_samvar_1_2.vcf.gz}_noE_samvar_1_2.vcf.NOsample.gz
        VCF_TUMORgermline=${SAMVCF%.vcf.gz}.GERMLINE.vcf
        VCF_TUMORsomatic=${SAMVCF%.vcf.gz}.SOMATIC.vcf
        if [ ! -f $VCF_TUMORgermline.gz ]; then
                command="java WriteDifference_TU_NO_asSomaticFilter $SAMVCF $NoVCF $VCF_TUMORgermline $VCF_TUMORsomatic; gzip $VCF_TUMORgermline ; gzip $VCF_TUMORsomatic"
                echo $command
                eval $command
        fi
done
# This is for both gatkUG and gatkHPCaller
for gatkVCF in `find $base_dir_wCNVs/alignments_vsn_k20/variants_default/ -name *gatk*.vcf | grep -v eval_1000 | grep -v GERMLINE | grep -v SOMATIC | grep -v combined | grep -v NO | grep -v paired `; do
        NoVCF=`echo ${gatkVCF} | sed 's/.gatk/.NOsample.gatk/'`
        VCF_TUMORgermline=${gatkVCF%.vcf}.GERMLINE.vcf
        VCF_TUMORsomatic=${gatkVCF%.vcf}.SOMATIC.vcf
        if [ ! -f $VCF_TUMORgermline ]; then
                command="java WriteDifference_TU_NO_asSomaticFilter $gatkVCF $NoVCF $VCF_TUMORgermline $VCF_TUMORsomatic"
                echo $command
                eval $command
        else
                echo "$VCF_TUMORgermline already exists"
        fi
done

# format HP Caller:
for vcf in `find $base_dir_wCNVs/alignments_*/variants_*/ -name *HPCaller*.vcf | grep -v rewritten | grep SOMATIC | grep -v paired`
do
        if [ ! -f ${vcf%.vcf}.rewritten.vcf ]; then
                command="java FormatHPCallerGATK_TUonly ${vcf} > ${vcf%.vcf}.rewritten.vcf"
                echo $command
                eval $command
        fi
done


# som sniper comma
for vcf in `find $base_dir_wCNVs/alignments_vsn_k20/variants_default/ -name *somaticSniper*.vcf | grep -v noComma | grep -v eval_1000 | grep -v GERMLINE | grep qual | grep -v combined | grep -v NO`; do
        outVCF=${vcf%.vcf}.noComma.vcf
        if [ ! -f $outVCF ]; then
                command="java Extract_Somatic_Mult_Alleles $vcf $outVCF"
                echo $command
                eval $command
        fi
done



