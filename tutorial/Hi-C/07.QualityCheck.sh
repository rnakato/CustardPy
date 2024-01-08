#!/bin/bash

build=hg38
gt=genometable.$build.txt

#sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.1.7.0.sif"
sing="singularity exec custardpy.sif"

outputdir=3DChromatin_ReplicateQC
mkdir -p $outputdir

samples="Control siCTCF siRad21" # Samples to be compared

chrs="chr21 chr22" # chromosomes to be considered
resolution=50000
norm=SCALE

# Generate metadatapairs
pairlist=$outputdir/metadata.pairs
rm -rf $pairlist
for sample1 in $samples; do for sample2 in $samples; do echo -e $sample1"\t"$sample2 >> $pairlist; done; done

# Generate contact data for all samples
rm -rf $outputdir/data
mkdir -p $outputdir/data
for cell in $samples
do
    hic=CustardPyResults_Hi-C/Juicer_hg38/$cell/aligned/inter_30.hic
    ### In case of starting from .hic files:
    #hic=hic/$cell/GSE196034_${cell}_merged.hic

    echo "preparing $cell..."
    for chr in $chrs; do
        $sing juicertools.sh dump observed $norm $hic $chr $chr BP $resolution \
        | awk -v chr=$chr 'OFS="\t" {printf("%s\t%d\t%s\t%d\t%d\n", chr, $1, chr, $2, $3)}' \
        | grep -v NaN > $outputdir/data/$cell.$chr.txt
    done

    cat $outputdir/data/$cell.*.txt > $outputdir/data/$cell.res$resolution
    $sing pigz $outputdir/data/$cell.res$resolution
    rm $outputdir/data/$cell.*.txt
done

# Generate samplelist
samplelist=$outputdir/metadata.samples
rm -rf $samplelist
for cell in $samples; do
    echo -e "$cell\t$(pwd)/$outputdir/data/$cell.res$resolution" >> $samplelist 
done

# Gnerate Binlist
binlist=$outputdir/data/Bins.$resolution.bed
rm -rf $binlist
for chr in $chrs; do
    $sing generate_binlist_from_gtfile.py $gt $chr $resolution >> $binlist
done
gzip -f $binlist

$sing run_3DChromatin_ReplicateQC.sh run_all \
      --metadata_samples $samplelist --bins $binlist.gz --metadata_pairs $pairlist --outdir $outputdir/output

# plot figure based on the results
$sing visualize_QC.py 3DChromatin_ReplicateQC/

### Optional
#$sing run_3DChromatin_ReplicateQC.sh concordance \
#      --metadata_samples $samplelist --bins $binlist.gz --metadata_pairs $pairlist --outdir $outputdir/output \
#      --methods HiCRep #GenomeDISCO

#$sing run_3DChromatin_ReplicateQC.sh summary \
#      --metadata_samples $samplelist --bins $binlist.gz --metadata_pairs $pairlist --outdir $outputdir/output \
#      --methods GenomeDISCO,HiCRep,HiC-Spector,QuASAR-Rep,QuASAR-QC

