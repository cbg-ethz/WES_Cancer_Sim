

fn_bcf=$1
organism=$2
fn_out=$3

if [ -z "$fn_bcf" ]; then 
	echo "usage $0 <fn_bcf> <organism> [fn_out]" 
	exit 0;
fi
if [ -z "$fn_out" ]; then 
	fn_out=${fn_bcf%.bcf}.sorted.bcf
fi

currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`
source `find ${gitDir} -name paths.sh`
source `find ${gitDir} -name genome.sh`

chrs=`$bcftools view $fn_bcf | grep -v "^#" | cut -f 1 | sort -u`

fn_tmp=/tmp/tmp.vcf

$bcftools view $fn_bcf | head -n 1000 | grep "^#" > $fn_tmp

for chr in $chrs; 
do
	n=`$bcftools view $fn_bcf | grep "^${chr}[[:space:]]" | wc -l`
	echo chr: $chr $n
	$bcftools view $fn_bcf | grep "^${chr}[[:space:]]" | sort -n -k2 >> $fn_tmp
done

nb=`$bcftools view $fn_bcf | wc -l`
na=`cat $fn_tmp | wc -l`

if [ ! "$na" -eq "$nb" ]; then 
	echo number of lines differ: before: $nb after: $na
	exit 1
fi

dict=`get_genome $organism`.dict

echo "$bcftools view -Sb $fn_tmp  -D $dict > $fn_out"
$bcftools view -Sb $fn_tmp  -D $dict > $fn_out
