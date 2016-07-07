
function get_genome()
{
	source `find $gitDir -name paths.sh`
	local organism=$1
	local ref=unknown
	if [ "$organim" == "9606" -o "$organism" == "human" ]; then
		#human genome
		ref=$humanRef
	fi
	echo $ref
}

function get_gtf()
{
	source `find $gitDir -name paths.sh`
	local organism=$1
	local ref=unknown
	if [ "$organism" == "9606" -o "$organism" == "human" ]; then
		#human genome
		ref=$humanGtf
	fi
	echo $ref
}


function get_gio()
{
	local organism=$1
	local ref=`get_genome $organism`
	local fn_gio=`dirname $ref`/genome.config
	if [ ! -f $fn_gio ]; then
		echo $make_gio $ref `dirname $ref` 1>&2
		$make_gio $ref `dirname $ref` 1>&2
	fi
	if [ ! -f $fn_gio ]; then
		echo Error: failed to compute the genome information object in $fn_gio
		exit 1; 
	fi
	echo $fn_gio
}
