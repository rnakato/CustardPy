#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname <build> <motifdir> <loop.bedpe>" 1>&2
}

if [ $# -ne 3 ]; then
  usage
  exit 1
fi

build=$1
motifdir=$2
loop=$3

ex(){ echo $1; eval $1; }

if test $build = "hg19" -o $build = "hg38" -o $build = "mm9" -o $build = "mm10"; then
    motiftext=/opt/motiffiles/$build.motifs.txt
else
    motiftext=""
fi

ex "java -Xms512m -Xmx32384m -jar /opt/juicer_tools.1.9.9_jcuda.0.8.jar motifs $build $motifdir $loop $motiftext"
