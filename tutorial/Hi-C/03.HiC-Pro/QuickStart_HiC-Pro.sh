sing="singularity exec --bind /work,/work2,/work3 --cleanenv /work/SingularityImages/custardpy.3.1.0.sif"

$sing gethg38genome.sh
genome=genome.hg38.fa

ncore=32

#sing="singularity exec custardpy.sif"

### Build Bowtie2 Index
indexdir=bowtie2-indexes
mkdir -p $indexdir
$sing bowtie2-build --threads $ncore $genome $indexdir/hg38
ln -rsf $genome $indexdir/hg38.fa

### Generate restriction site file
resfile=MboI_resfrag_hg38.bed
$sing /usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py -r mboi -o $resfile $genome

### Edit the config file (config-hicpro.txt)
# Please edit the lines BOWTIE2_IDX_PATH, GENOME_SIZE, and GENOME_FRAGMENT in config-hicpro.txt to match your directory structure.
# Note that you need to specify absolute paths for these entries.

### Run Hi-C-Pro
configfile=config-hicpro.txt
hicprodir=HiCProResults

$sing hicpro -i ../fastq/ -o $hicprodir -c $configfile
