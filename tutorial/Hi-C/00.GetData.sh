#mkdir -p fastq/Hap1-A
#fastq-dump --split-files --gzip SRR5266584 -O fastq/Hap1-A

# insteadn, download the paired fastq files from our Google Drive:
# https://drive.google.com/file/d/1Clo_sh3n4r82lp3ynkmUdLzcGx3gpu5Q/view?usp=sharing


sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.1.0.sif"

$sing gethg38genome.sh

genome=genome.hg38.fa
indexdir=bwa-indexes
mkdir -p $indexdir
$sing bwa index -p $indexdir/hg38 $genome
ln -rsf $genome $indexdir/hg38
