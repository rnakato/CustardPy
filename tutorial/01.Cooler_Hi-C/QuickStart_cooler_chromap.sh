#!/bin/bash

build=hg38
gt=genometable.$build.txt
index_chromap=chromap-indexes/$build
genome=genome.$build.fa
ncore=64
enzyme=MboI
gene=refFlat.$build.txt
fastq_post="_"  # "_" or "_R"
qthre=10

sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.5.2.sif"
sing_gpu="apptainer exec --nv --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.5.2.sif"

for cell in Control siCTCF siRad21 siNIPBL
do
    # generate .cool and .hic files
    $sing custardpy_cooler -M chromap -o CustardPyResults_chromap \
          -g $gt -f $genome -b $build -e $enzyme -z $fastq_post -q $qthre \
          -i $index_chromap -p $ncore fastq/$cell $cell

    # downstream analysis of .cool files (eigenvector, insulation score, loop calling, etc.)
    odir=CustardPyResults_chromap/Cooler_$build/$cell
    cool=$odir/cool/contact.chromap.q30.mcool
    $sing custardpy_process_cool -t $ncore -g $gt -f $genome $cool $odir $cell

    ## (Optional) loop calling with fithic (take long time, so run separately)
    resolutions_fithic=25000
    odir_fithic=$odir/loops/fithic/$resolutions_fithic
    $sing run_fithic.sh -g $gt $cool $odir_fithic $cell $resolutions_fithic \
    | tee $odir/log/fithic.$resolutions_fithic.log

    ## (Optional) Juicer analysis (requires GPU, so run separately)
    hic=$odir/hic/contact.chromap.q30.hic
    norm=SCALE
    $sing_gpu custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
done
