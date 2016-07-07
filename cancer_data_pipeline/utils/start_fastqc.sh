#/bin/bash -l

dir_=$dir_
tsv_dir=$tsv_dir
out=$out

echo "$dir_"
echo "$tsv_dir"
echo "$out"

if [ -z $dir_ ]; then
	echo "Could not find the directory dir_."
	exit -1
fi
if [ -z $tsv_dir ]; then
        echo "Could not find the directory of the tsv-file."
        exit -1
fi
if [ -z $out ]; then
        echo "Could not find the output directory."
        exit -1
fi


currScriptDir=`dirname $0`
source `find ${currScriptDir}/../../ -name paths.sh`
source `find ${dir_} -name utils.sh`
source `find ${dir_} -name paths.sh`

run_fastqc $tsv_dir $out $NSLOTS
