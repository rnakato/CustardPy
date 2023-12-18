# Changelog

## CustardPy

### 1.5.0 (2023-12-18)
- AluI restriction enzyme added to Juicer and Cooler

### 1.4.5 (2023-12-01)
- Removed LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/compat/:/usr/local/cuda/lib64
- Bug fix for HICCUPS error (CUDA not detected)

### 1.4.4 (2023-11-28)
- Added `-o` option to `custardpy_process_hic`, `custardpy_cooler_subfunc.sh`, `makeMatrix_intra.sh`, `call_HiCCUPS.sh` and `juicer_callTAD.sh`

### 1.4.3 (2023-11-05)
- Added Juicerstats.sh to summarize the mapping stats of Juicer
- Added Coolerstats.sh to summarize the mapping stats of Cooler

### 1.4.2 (2023-10-24)
- Added LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/compat/:/usr/local/cuda/lib64

### 1.4.1 (2023-9-23)
- Bug fix in installaion of 3DChromatin_ReplicateQC and HiCrep
- Bug fix in QuASAR (https://github.com/kundajelab/3DChromatin_ReplicateQC/issues/13)

### 1.4.0 (2023-8-01)
- Bug fix for HOMER installation in Dockerfile
- Changed the parameter of BWA MEM from `-SP` to `-5SP`

### 1.3.1 (2023-7-08)
- Changed the name of output directory for `custardpy_juicer`, `custardpy_cooler_HiC`, and `custardpy_cooler_MicroC`

### 1.3.0 (2023-7-06)
- Modified `custardpy_cooler_HiC` and `custardpy_cooler_MicroC` to allow multiple FASTQ files
- Changed the name of output directory for `custardpy_juicer`, `custardpy_cooler_HiC`, and `custardpy_cooler_MicroC`
- Modified options of Pairtools according to the [tutorial](https://github.com/open2c/pairtools/blob/master/doc/examples/pairtools_walkthrough.ipynb)
- Added [python-FIREcaller](https://github.com/jakublipinski/python-FIREcaller)
- Added the restriction files of Cooler for various genome builds and enzymes (see the help of `custardpy_cooler_HiC`)
- Bug fix: Wrong PATH of pairsqc.py

### 1.2.0 (2023-6-11)
- Added [3DChromatin_ReplicateQC](https://github.com/kundajelab/3DChromatin_ReplicateQC)
- Added `run_3DChromatin_ReplicateQC.sh` for executing 3DChromatin_ReplicateQC
- Updated `custardpy_cooler_HiC` and `custardpy_cooler_MicroC`

### 1.1.1 (2023-06-07)
- Added [GENOVA](https://github.com/robinweide/GENOVA)
- Added `calculate_compartment_strength` command to calculate compartment strength
- Added [HiCUP](https://www.bioinformatics.babraham.ac.uk/projects/hicup/read_the_docs/html/)
- Bug fix of pairtools and numpy version conflict

### 1.1.0 (2023-05-24)
- Added `DEG_boundary_analysis` command for TAD-proximity analysis
- Added `custardpy_clustering_boundary` command for clustering TAD boundaries
- Added `custardpy_differential_DRF` command for differential DRF analysis
- Added `calculate_compartment_strength` command for calculating compartment strength
- Bug fix in `makeEigen.sh`
- Changed name of `custardpy_mappingHiC` to `custardpy_cooler_HiC`
- Changed name of `custardpy_mappingMicroC` to `custardpy_cooler_MicroC`

### 1.0.0 (2023-05-17)
- Major Release!
- Merged the GitHub repositories for CustardPy (this site) and [Docker_CustardPy](https://github.com/rnakato/Docker_CustardPy) into a single repository.
- Unified the Docker images for [CustardPy](https://hub.docker.com/r/rnakato/custardpy) and [CustardPy_Juicer](https://hub.docker.com/r/rnakato/custardpy_juicer). Version 1 of the [CustardPy](https://hub.docker.com/r/rnakato/custardpy) docker image now supports all analyses previously offered by CustardPy and CustardPy_Juicer, rendering the latter unnecessary.

- Docker Image Specifications:
    - OS: Ubuntu 20.04
    - CUDA: cuda11.0.3
    - cuDNN: v8
    - R: 4.3
    - Python: 3.9.13
    - OpenJDK v17
    - Juicer v1.6
    - Juicer tools v1.22.01

Note: CustardPy does not support Juicer Tools v2 due to incompatibility with the ``.cool`` file format.

- Added tools
    - [GENOVA](https://github.com/robinweide/GENOVA)
    - [CHESS](https://chess-hic.readthedocs.io/en/latest/index.html)
    - [FIREcaller](https://github.com/yycunc/FIREcaller)
    - [STRIPENN](https://github.com/VahediLab/stripenn)
    - [HiC-Pro](https://github.com/nservant/HiC-Pro)
    - [coolpup.py](https://github.com/open2c/coolpuppy)
    - [FitHiChIP](https://ay-lab.github.io/FitHiChIP/html/index.html)
    - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

- Merged `drawTriangleRatioMulti` to `plotHiCfeature`
- Fixed the error "QStandardPaths: XDG_RUNTIME_DIR points to non-existing path '/run/user/1000', please create it with 0700 permissions." when calculating PC1.
- Removed /root/.cpanm/work directory

### 0.4.1 (2023-02-24)
- Added `libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-render-util0` libraries to avoid an error 'qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.'

### 0.4.0 (2023-01-26)
- Updated R from 3.6 to 4.2
- Added [FitHiC](https://github.com/ay-lab/fithic)
- Added [ChIAPop](https://github.com/wh90999/ChIAPoP) for ChIA-PET analysis

### 0.3.0 (2022-10-28)
- Downgraded pairtools from v1.0.1 to v0.3.0

<!--
- change custardpy_mappingMicroC not to output BAM file that takes long time
- Added `--backend cython` option to *pairtools dedup* for the consistensy to pairtools v0.3.0
-->

### 0.2.0 (2022-08-31)
- Public release
- Bug Fix for visualization
- Updated Manual


### 0.1.0
- First commit



## Docker_CustardPy_Juicer (deprecated)

### 0.4.0 (2023-04-25)
- Changed OS from Ubuntu16.06 to Ubuntu20.04
- Changed CUDA driver from v8 to 11.0.3

### 0.3.2 (2023-02-07)
- Change the upper limit of memory in Java from `-Xmx32768m` to `-Xmx128768m`

### 0.3.1 (2022-12-12)
- Modify mega.sh to adjust CustardPy_juicer
- Added juicer_genhic.sh that generates a .hic file from the mapped data
- Deprecated `juicer_postprocessing.sh`

### 0.3.0 (2022-10-27)
- Downgrade `juicer_tools` from juicer_tools.2.13.07.jar to juicer_tools.1.22.01.jar to keep the consistency between .hic and .cool (see https://github.com/deeptools/HiCExplorer/issues/798)

### 0.2.0 (2022-08-31)
- Added restriction sites for Arima and Sau3AI
- Added restriction sites for various species
- change R installation from conda to apt

### 0.1.0
- First commit
