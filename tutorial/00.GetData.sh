mkdir -p fastq/Hap1-A
fastq-dump --split-files --gzip SRR5266584 -O fastq/Hap1-A

# insteadn, download the paired fastq files from our Google Drive:
# https://drive.google.com/file/d/1Clo_sh3n4r82lp3ynkmUdLzcGx3gpu5Q/view?usp=sharing
