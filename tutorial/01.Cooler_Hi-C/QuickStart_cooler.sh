#!/bin/bash

build=hg38
gt=genometable.$build.txt
index_bwa=bwa-indexes/$build
genome=genome.$build.fa
ncore=64
enzyme=MboI
gene=refFlat.$build.txt
fastq_post="_"  # "_" or "_R"

version=3.5.1
sing="apptainer exec --nv --bind /work,/work2,/work3 /work/SingularityImages/custardpy.$version.sif"

for cell in Control siCTCF siRad21 siNIPBL
do
    # generate .cool and .hic files
    $sing custardpy_cooler -o CustardPyResults_$version -g $gt -f $genome -b $build -e $enzyme -z $fastq_post -i $index_bwa -p $ncore fastq/$cell $cell

    # downstream analysis of .cool files (eigenvector, insulation score, loop calling, etc.)
    odir=CustardPyResults_$version/Cooler_$build/$cell
    cool=$odir/cool/contact.bwa.q30.mcool
    $sing custardpy_process_cool -t $ncore -g $gt -f $genome $cool $odir $cell

    ## (Optional) loop calling with fithic (take long time, so run separately)
    resolutions_fithic=25000
    odir_fithic=$odir/loops/fithic/$resolutions_fithic
    $sing run_fithic.sh -g $gt $cool $odir_fithic $cell $resolutions_fithic \
       | tee $odir/log/fithic.$resolutions_fithic.log

    mkdir -p $odir/loops/mustache
    $sing mustache -f $cool -r 5kb -norm weight -pt 0.05 -p $ncore \
            -o $odir/loops/mustache/loop.cool.5kb.tsv


    ## (Optional) Juicer analysis (requires GPU, so run separately)
    hic=$odir/hic/contact.bwa.q30.hic
    norm=SCALE
    $sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
done
