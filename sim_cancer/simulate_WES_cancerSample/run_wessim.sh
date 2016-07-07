currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`
source `find ${gitDir} -name paths.sh`

if [ $# -lt 3 ]; then
        echo "found $# arguments"
        echo "Usage $0 <fn_fasta> <fn_probes> <out_dir> <num_reads>"
        echo "options: "
        echo "    --read-len <int>    :(default 100)"
        echo "    --skip-after-blat   :only do primary data processing (default no)"
        echo "    --num-cpus              :(default 1)"
        exit 0
fi

# arguments
fn_ref=$1; shift
fn_probes=$1; shift
out_dir=$1; shift
num_reads=$1; shift
read_len=100
skip_after_blat=no
num_cpus=1
# defaults for the blat alignments of probes against genomes
blat_min_identity=90 # same as the wessim default
blat_min_score=100 # same as the wessim default
blat_port=6666

while [ 0 -lt $# ]; do
    if [ "$1" == "--read-len" ]; then
                shift
        read_len=$1
                echo use read length: $read_len
        elif [ "$1" == "--skip-after-blat" ]; then
                skip_after_blat=yes
        elif [ "$1" == "--num-cpus" ]; then
                shift
                num_cpus=$1
        elif [ "$1" == "--blat-min-identity" ]; then
                shift
                blat_min_identity=$1
        elif [ "$1" == "--blat-min-score" ]; then
                shift
                blat_min_score=$1
        elif [ "$1" == "--blat-port" ]; then
                shift
                blat_port=$1
        else
                echo did not understand $1
                echo will be ignored...
    fi
    shift
done

#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/gfClient
#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/gfServer
#http://sourceforge.net/projects/gemsim/files/latest/download (by hand) 

tools_dir=$gitDir/sim_cancer/tools/
faToTwoBit=$tools_dir/faToTwoBit 
gfServer=$tools_dir/gfServer
gfClient=$tools_dir/gfClient
gemsim_models=$tools_dir/GemSIM_v1.6/models/
fn_model=$tools_dir/GemSIM_v1.6/models/ill100v5_p.gzip
wessim_dir=$tools_dir/Wessim_beta

mkdir -p $out_dir
if [ ! -f $out_dir/`basename $fn_ref` ]; then 
	ln -s $fn_ref $out_dir
fi
fn_ref=$out_dir/`basename $fn_ref`
fn_blat=$out_dir/probe.txt.fa.psl

if [ ! -f $fn_ref.fai ]; then 
	echo $samtools faidx $fn_ref
	$samtools faidx $fn_ref
fi
if [ ! -f ${fn_ref%fa}2bit ]; then 
	echo $faToTwoBit $fn_ref ${fn_ref%fa}2bit
	$faToTwoBit $fn_ref ${fn_ref%fa}2bit
fi


# Generate a FASTA file of probe sequence
if [ ! -f $fn_probes.fa ]; then 
	${localToolsDir}/pipeline/bin/python2.7 $wessim_dir/Prep_Probe2Fa.py $fn_probes
fi

if [ ! -f $fn_blat ]; then 
	# Establish your local blat server
	#server_running=`ps aux | grep gfServer | grep -v grep`
	echo "$gfServer start -canStop localhost $blat_port ${fn_ref%fa}2bit &"
	$gfServer start -canStop localhost $blat_port ${fn_ref%fa}2bit & #(You need to provide the abosolute path of the reference file. gfServer will consume one whole thread, you need to use a separated thread to continue the following steps)
	# wait for blat server to read reference
	echo wait for blat server ...
	sleep 300
	# Run blat search to generate the match list
	echo "PATH=`dirname $gfClient`:$PATH ${localToolsDir}/pipeline/bin/python2.7 $wessim_dir/Prep_BlatSearch.py -R ${fn_ref%fa}2bit -P $fn_probes.fa -p $blat_port -o $fn_blat"
	PATH=`dirname $gfClient`:$PATH ${localToolsDir}/pipeline/bin/python2.7 $wessim_dir/Prep_BlatSearch.py -R ${fn_ref%fa}2bit -P $fn_probes.fa -p $blat_port -o $fn_blat -i $blat_min_identity -s $blat_min_score #(Note that the path to ref.2bit is not based on your local machine. It should be used without path, because the gfServer has it in its root)

	# stop the blat server
	$gfServer stop localhost $blat_port
fi
# Run Wessim2 in probe hybridization mode.
# This will generate *result_1.fastq.gz* and *result_2.fastq.gz* (paired-end mode / gzip compressed).
if [ ! "$skip_after_blat" == "yes" ]; then 
	if [ ! -f $out_dir/result_1.fastq.gz ]; then 
		cd $wessim_dir
		${localToolsDir}/pipeline/bin/python2.7 Wessim2.py -R $fn_ref -t $num_cpus -P $fn_probes.fa -B $fn_blat -n $num_reads -l $read_len -M $fn_model -pz -o $out_dir/result
	fi
fi


