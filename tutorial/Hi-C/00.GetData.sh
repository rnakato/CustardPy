### Download FASTQ files
mkdir -p fastq/Hap1-A fastq/WaplKO_3.3-A
fastq-dump --split-files --gzip SRR5266584 -O fastq/Hap1-A
fastq-dump --split-files --gzip SRR5266585 -O fastq/WaplKO_3.3-A

# instead, download the paired fastq files from our Google Drive:
# https://drive.google.com/file/d/1GA3d-TD7RpYFhP4VYcnN1WQmV3FjC58N/view?usp=sharing
# https://drive.google.com/file/d/1GDB1Y35heroYOHg1CGooQc8BF25jKZYc/view?usp=sharing


#sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.2.0.sif"
sing="singularity exec custardpy.sif"

$sing gethg38genome.sh

genome=genome.hg38.fa
indexdir=bwa-indexes
mkdir -p $indexdir
$sing bwa index -p $indexdir/hg38 $genome
ln -rsf $genome $indexdir/hg38

indexdir=chromap-indexes
mkdir -p $indexdir
chromap -i -t 12 -r $genome -o $indexdir/hg38
