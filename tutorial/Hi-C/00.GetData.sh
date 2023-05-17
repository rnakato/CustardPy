### Download FASTQ files
mkdir -p fastq/Hap1-A fastq/WaplKO_3.3-A
fastq-dump --split-files --gzip SRR5266584 -O fastq/Hap1-A
fastq-dump --split-files --gzip SRR5266585 -O fastq/WaplKO_3.3-A

# Instead, download the paired fastq files from our website:
# wget -nv --timestamping https://onl.la/CJi7U3m -O Hap1-A.zip
# wget -nv --timestamping https://onl.la/93jAHRe -O WaplKO_3.3-A.zip


#sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.1.0.0.sif"
sing="singularity exec custardpy.sif"

$sing gethg38genome.sh

genome=genome.hg38.fa
indexdir=bwa-indexes
mkdir -p $indexdir
$sing bwa index -p $indexdir/hg38 $genome
ln -rsf $genome $indexdir/hg38

indexdir=chromap-indexes
mkdir -p $indexdir
$sing chromap -i -t 12 -r $genome -o $indexdir/hg38
