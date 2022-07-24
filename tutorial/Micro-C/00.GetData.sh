#mkdir -p fastq
#fastq-dump --split-files --gzip SRR8954797 -O fastq/

# instead, download the paired fastq files from our Google Drive:
# https://drive.google.com/drive/folders/16gn5uO1u9yc5wqh7XQHg844FiYqzrnTB?usp=sharing

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.1.0.sif"

$sing getmm39genome.sh

genome=genome.mm39.fa
indexdir=bwa-indexes
mkdir -p $indexdir
$sing bwa index -p $indexdir/mm39 $genome
ln -rsf $genome $indexdir/mm39.fa
