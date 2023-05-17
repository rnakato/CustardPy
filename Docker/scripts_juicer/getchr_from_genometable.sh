#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname <genome table>" 1>&2
}

if [ $# -ne 1 ]; then
  usage
  exit 1
fi

gt=$1
cut -f1 $gt | tr '\n' ' '
