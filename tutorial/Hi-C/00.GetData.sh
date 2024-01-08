### In case of starting from FASTQ files: Download FASTQ files
mkdir -p fastq/siCTCF fastq/siRad21 fastq/Control # fastq/siNIPBL

# siCTCF
wget -nv --timestamping ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/013/SRR17870713/SRR17870713_1.fastq.gz -P fastq/siCTCF
wget -nv --timestamping ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/013/SRR17870713/SRR17870713_2.fastq.gz -P fastq/siCTCF
# Control
wget -nv --timestamping ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/018/SRR17870718/SRR17870718_1.fastq.gz -P fastq/Control
wget -nv --timestamping ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/018/SRR17870718/SRR17870718_2.fastq.gz -P fastq/Control
# siRad21
wget -nv --timestamping ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/040/SRR17870740/SRR17870740_1.fastq.gz -P fastq/siRad21
wget -nv --timestamping ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/040/SRR17870740/SRR17870740_2.fastq.gz -P fastq/siRad21
# siNIPBL
#wget -nv --timestamping ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/029/SRR17870729/SRR17870729_1.fastq.gz -P fastq/siNIPBL
#wget -nv --timestamping ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/029/SRR17870729/SRR17870729_2.fastq.gz -P fastq/siNIPBL

### In case of starting from .hic files: download .hic files
mkdir -p hic/siCTCF hic/siRad21 hic/Control # hic/siNIPBL
wget -nv --timestamping https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196034/suppl/GSE196034%5FCTCFKD%5Fmerged.hic -P hic/siCTCF
wget -nv --timestamping https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196034/suppl/GSE196034%5FRad21KD%5Fmerged.hic -P hic/siRad21
wget -nv --timestamping https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196034/suppl/GSE196034%5FControl%5Fmerged.hic -P hic/Control
#wget -nv --timestamping https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196034/suppl/GSE196034%5FNIPBLKD%5Fmerged.hic -P hic/siNIPBL


### Download the reference genome and build genome index
#sing="singularity exec --nv --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.1.5.0.sif"
sing="singularity exec custardpy.sif"

$sing gethg38genome.sh

genome=genome.hg38.fa
indexdir=bwa-indexes
mkdir -p $indexdir
$sing bwa index -p $indexdir/hg38 $genome
ln -rsf $genome $indexdir/hg38

#indexdir=chromap-indexes
#mkdir -p $indexdir
#$sing chromap -i -t 12 -r $genome -o $indexdir/hg38
