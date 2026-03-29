#!/bin/bash

build=hg38
gt=genometable.$build.txt
index_bwa=bwa-indexes/$build
genome=genome.$build.fa
ncore=64
enzyme=MboI
fastq_post="_"  # "_" or "_R"

sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.2.0.sif"
sing_gpu="apptainer exec --nv --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.2.0.sif"

for cell in Control siCTCF siRad21 siNIPBL
do
    # generate .cool and .hic files
    $sing custardpy_cooler -g $gt -f $genome -b $build -e $enzyme -z $fastq_post -i $index_bwa -p $ncore fastq/$cell $cell

    # downstream analysis of .cool files (eigenvector, insulation score, loop calling, etc.)
    odir=CustardPyResults/Cooler_$build/$cell
    cool=$odir/cool/cooler.dedup.q30.multires.cool
    $sing custardpy_process_cool -t $ncore -g $gt -f $genome $cool $odir $cell
#exit
    ## (Optional) loop calling with fithic (take long time, so run separately)
    resolutions_fithic=25000
    pairfile=$odir/pairs/dedup.bwa.q30.pairs.gz
    odir_fithic=$odir/loops/fithic/$resolutions_fithic
    $sing run_fithic.sh -g $gt $pairfile $odir_fithic $cell $resolutions_fithic | tee $odir/log/fithic.$resolutions_fithic.log

    ## (Optional) Juicer analysis (requires GPU, so run separately)
    hic=$odir/hic/contact_map.q30.hic
    norm=SCALE
#    $sing_gpu custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
done
