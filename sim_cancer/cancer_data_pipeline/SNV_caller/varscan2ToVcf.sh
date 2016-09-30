#/bin/bash -l

source `find ${gitDir} -name paths.sh`
source `find ${gitDir} -name genome.sh`
source `find ${gitDir} -name utils.sh`


# DESCRIPTION_______________________________________________________________________________________________________________________________________________________________________________
# This script converts the snv list into the vcf format.
#___________________________________________________________________________________________________________________________________________________________________________________________

if [ -z "$1" -o -z "$2" ]; then
	echo "Usage: $0 <SNVs_varscan2.txt.snp.Somatic.hc> <fn_genome>"
	exit -1
fi


varscan_list=$1
fn_genome=$2
varscanVCF=${varscan_list}.vcf


echo -e "##fileformat=VCFv4.0" > $varscanVCF
echo -e "##source=VarScan2" >> $varscanVCF
echo -e "##reference=$fn_genome" >> $varscanVCF
echo -e "##INFO=<ID=NR,Number=1,Type=Integer,Description="Number of reads supporting reference in normal">" >> $varscanVCF
echo -e "##INFO=<ID=NV,Number=1,Type=Integer,Description="Number of reads supporting variant in normal">" >> $varscanVCF
echo -e "##INFO=<ID=TR,Number=1,Type=Integer,Description="Number of reads supporting reference in tumor">" >> $varscanVCF
echo -e "##INFO=<ID=TV,Number=1,Type=Integer,Description="Number of reads supporting variant in tumor">" >> $varscanVCF
echo -e "##INFO=<ID=NVF,Number=.,Type=Float,Description="Variant frequency in normal sample \(in \%\)">" >> $varscanVCF
echo -e "##INFO=<ID=TVF,Number=.,Type=Float,Description="Variant frequency in tumor \(in \%\)">" >> $varscanVCF
echo -e "##INFO=<ID=NG,Number=1,Type=String,Description="Normal genotype">" >> $varscanVCF
echo -e "##INFO=<ID=TG,Number=1,Type=String,Description="Tumor genotype">" >> $varscanVCF
echo -e "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="If status of variant is somatic">" >> $varscanVCF
echo -e "##INFO=<ID=VPV,Number=.,Type=Float,Description="Variant p-value: Significance of variant read count vs. baseline error rate">" >> $varscanVCF
echo -e "##INFO=<ID=SPV,Number=.,Type=Float,Description="Somatic p-value: Significance of tumor read count vs. normal read count">" >> $varscanVCF
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> $varscanVCF


#__ Check how many SNVs were found ________________________________________________________________________________________________________________________________________________________
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
nr_snvs=$(cat $varscan_list | awk 'NR > 1' | wc -l)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#__ Loop over all of them and write them into the vcf format ______________________________________________________________________________________________________________________________

first=true
while read line ; 
do
	if $first; then 
		first=false
		continue
	fi	

	#__ Extract chromosome, position, ... _____________________________________________________________________________________________________________________________________________
	chr=$( echo $line | awk '{print $1}') 
	pos=$( echo $line | awk '{print $2}')
	ref=$( echo $line | awk '{print $3}')
	var=$( echo $line | awk '{print $4}')
	normal_reads1=$( echo $line | awk '{print $5}')
	normal_reads2=$( echo $line | awk '{print $6}')
	normal_var_freq=$( echo $line | awk '{print $7}')
	normal_gt=$( echo $line | awk '{print $8}')
	tumor_reads1=$( echo $line | awk '{print $9}')
	tumor_reads2=$( echo $line | awk '{print $10}')
	tumor_var_freq=$( echo $line | awk '{print $11}')
	tumor_gt=$( echo $line | awk '{print $12}')
	somatic_status=$( echo $line | awk '{print $13}')
	variant_p_value=$( echo $line | awk '{print $14}')
	somatic_p_value=$( echo $line | awk '{print $15}')

	normal_var_freq_percent=${normal_var_freq%\%}
	tumor_var_freq_percent=${tumor_var_freq%\%}

	if [ $somatic_status == "Somatic" ]; then
		INFO=$(echo "NR=$normal_reads1;NV=$normal_reads2;TR=$tumor_reads1;TV=$tumor_reads2;NVF=$normal_var_freq_percent;TVF=$tumor_var_freq_percent;NG=$normal_gt;TG=$tumor_gt;SOMATIC;VPV=$variant_p_value;SPV=$somatic_p_value")
	else
		INFO=$(echo "NR=$normal_reads1;NV=$normal_reads2;TR=$tumor_reads1;TV=$tumor_reads2;NVF=$normal_var_freq_percent;TVF=$tumor_var_freq_percent;NG=$normal_gt;TG=$tumor_gt;VPV=$variant_p_value;SPV=$somatic_p_value")
	fi

	echo -e "$chr\t$pos\t.\t$ref\t$var\t.\t.\t$INFO" >> $varscanVCF

done < $varscan_list
#___________________________________________________________________________________________________________________________________________________________________________________________
