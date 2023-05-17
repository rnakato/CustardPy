mkdir -p fastq
fastq-dump --split-files --gzip SRR8954797 -O fastq/

# Instead, download the paired fastq files from our Google Drive:
# wget -nv --timestamping hhttps://onl.la/6wy4EKa -O fastq/SRR8954797_1.fastq.gz
# wget -nv --timestamping https://onl.la/QbvvNLB -O fastq/SRR8954797_2.fastq.gz

#sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.1.0.0.sif"
sing="singularity exec custardpy.sif"

$sing getmm39genome.sh

genome=genome.mm39.fa
indexdir=bwa-indexes
mkdir -p $indexdir
$sing bwa index -p $indexdir/mm39 $genome
ln -rsf $genome $indexdir/mm39

indexdir=chromap-indexes
mkdir -p $indexdir
$sing chromap -i -t 12 -r $genome -o $indexdir/mm39
