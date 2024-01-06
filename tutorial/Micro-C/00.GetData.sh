mkdir -p fastq/C36_rep1

### In case of downloading FASTQ files using fastq-dump
fastq-dump --split-files --gzip SRR16763198 -O fastq/C36_rep1

### In case of downloading FASTQ files from ENA (https://www.ebi.ac.uk/ena/)
wget -nv --timestamping ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/098/SRR16763198/SRR16763198_1.fastq.gz -P fastq/C36_rep1
wget -nv --timestamping ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/098/SRR16763198/SRR16763198_2.fastq.gz -P fastq/C36_rep1

### Download the reference genome and build genome index
#sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.1.6.0.sif"
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
