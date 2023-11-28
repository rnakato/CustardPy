#!/bin/bash

cmdname=`basename $0`
function usage()
{
    echo "$cmdname [-o] <command>" 1>&2
    echo '   <command>: Command of juicer_tools.jar' 1>&2
    echo '   Options:' 1>&2
    echo '     -o: Use older version of juicer_tools.jar (juicer_tools.1.9.9_jcuda.0.8.jar, default: juicer_tools.1.22.01.jar)' 1>&2
}

useoldversion="no"
while getopts o option; do
    case ${option} in
        o) useoldversion="yes" ;;
        \?) 
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        *)
            usage
            exit 1
            ;;
    esac
done
shift $((OPTIND - 1))

if [ $# -eq 0 ]; then
  usage
  exit 1
fi

if test $useoldversion = "no"; then
    # juicer_tools.1.22.01.jar
    command="java -Xms512m -Xmx64384m -jar /opt/juicer/scripts/common/juicer_tools.jar $@"
else 
    command="java -Xms512m -Xmx64384m -jar /opt/juicer/juicer_tools.1.9.9_jcuda.0.8.jar $@"
fi

echo $command
eval $command
