#!/bin/bash

cmdname=`basename $0`
function usage()
{
    echo "$cmdname: Run command for paired fastq files in all directories within the fastq directory."
    echo "Usage: $cmdname [-p [1|2]] <fastq directory> <command>"
    echo "  <fastq directory> Directory containing paired FASTQ files (e.g., 'fastq/')"
    echo "  <command>  specify the command to execute (should be in quotes)"
    echo "  -p postfix"
    echo "       1: '*_1.fastq.gz' and '*_2.fastq.gz' (default)"
    echo "       2: '*_R1.fastq.gz' and '*_R2.fastq.gz'"
    echo "  -s delimiter"
    echo "       s: use space (default)"
    echo "       c: use comma"
}

postfix1=_1.fastq.gz
postfix2=_2.fastq.gz
delimiter=" "
while getopts d:p:s: option; do
    case ${option} in
        d) fqdir=${OPTARG} ;;
        p) case ${OPTARG} in
            1)
                postfix1=_1.fastq.gz
                postfix2=_2.fastq.gz
                ;;
            2)
                postfix1=_R1.fastq.gz
                postfix2=_R2.fastq.gz
                ;;
            *)
                echo "Error: Specify 1 or 2 to the option '-p'."
                usage
                exit 1
                ;;
           esac
           ;;
        s) case ${OPTARG} in
            s) delimiter=" " ;;
            c) delimiter="," ;;
            *)
                echo "Error: Specify 's' or 'c' to the option '-s'."
                usage
                exit 1
                ;;
           esac
           ;;
        \?) 
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        *)
        echo "test"
            usage
            exit 1
            ;;
    esac
done
shift $((OPTIND - 1))

if [ $# -ne 2 ]; then
    usage
    exit 0
fi

fqdir=$1
command="$2"

if [ ! -d "$fqdir" ]; then
    echo "Error: $fqdir is not a directory."
    exit
fi
for dir in $fqdir/*; do
    if [ ! -d "$dir" ]; then continue; fi
    
    fq1_list=()
    fq2_list=()

    for fq1 in "$dir"/*"$postfix1"; do
        if [ ! -f "$fq1" ]; then
            echo "$dir: *${postfix} does not exist. Skipping"
            continue
        fi

        fq2="${fq1%$postfix1}$postfix2"
        if [ ! -f "$fq2" ]; then 
            echo "$dir: $fq2 does not exist. Skipping"
            continue
        fi

        fq1_list+=("$fq1")
        fq2_list+=("$fq2")
    done

    fq1_list=$(IFS=$delimiter; echo "${fq1_list[*]}")
    fq2_list=$(IFS=$delimiter; echo "${fq2_list[*]}")

    if [ -n "$fq1_list" ] && [ -n "$fq2_list" ]; then
        $command $fq1_list $fq2_list
    fi
done
