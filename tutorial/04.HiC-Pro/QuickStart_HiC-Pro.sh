sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.1.0.sif"

### Edit the config file (config-hicpro.txt)
# Please edit the lines BOWTIE2_IDX_PATH, GENOME_SIZE, and GENOME_FRAGMENT in config-hicpro.txt to match your directory structure.
# Note that you need to specify absolute paths for these entries.

### Run Hi-C-Pro
configfile=config-hicpro.txt
hicprodir=HiCProResults

$sing hicpro -i fastq/ -o $hicprodir -c $configfile
