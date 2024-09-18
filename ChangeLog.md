# Changelog

## CustardPy

### 2.1.2 (2024-9-18)
- Fixed a bug where `run_3DChromatin_ReplicateQC.sh` gave a "No such file or directory" error.

### 2.1.1 (2024-9-12)
- Modified `convert_JuicerDump_to_dense.py` in **CustardPy** to fix a future warning message.

### 2.1.0 (2024-9-10)
- Fixed a bug where ``juicer_map.sh`` failed when passing "\_" for ``fastq_post`` valuable and there were underscores (\_) in the FASTQ directory name.

### 2.0.1 (2024-8-27)
- Added `libgit2-dev` to install `gert` in R
- Added a message in ``juicer_pigz.sh`` and ``juicer_unpigz.sh``

### 2.0.0 (2024-5-22)
- Changed Python environment from conda to micromamba (`/opt/micromamba`)
- Updated Python from 3.8.17 to 3.10.13
- Install [HiCLift](https://github.com/XiaoTaoWang/HiCLift)

### 1.9.1 (2024-04-04)
- Updated  HiC1Dmetrics from v0.2.9 to v0.2.10
- Omitted `${juiceDir}/scripts/common/check.sh` which checks the number of ${outputdir}/merged_sort.txt, ${outputdir}/merged_nodups.txt, ${outputdir}/dups.txt and ${outputdir}/opt_dups.txt.

### 1.9.0 (2024-03-21)
- Updated `distance_vs_count.Juicer` and `distance_vs_count.Juicer.log` to use minimum Q-value threshold of 30 and window size 50 kbp by default
- Updated `plot_distance_count.R` to set the background color to white
- Updated `plot_distance_count.log.R` to set the background color to white and limit the x- and y-axes
- Added `plot_distance_count_all.R`, `plot_distance_count_all.log.R`, `plot_distance_count_multi.R`, `plot_distance_count_multi.log.R`
- Modified `juicer_pigz.sh` and `juicer_unpigz.sh` to compress/uncompress aligned/dups.txt, aligned/opt_dups.txt and aligned/collisions*txt as well
- `makeEigen.sh`: Fixed the bug where the computation would not finish if the computation of PC1 in h1d failed.
- Fixed a bug where `pairtools` was not found.
- Added [CALDER2](https://github.com/CSOgroup/CALDER2)
- Updated chromap from v0.2.4 to v0.2.6
- Updated bedtools from v2.30.0 to v2.31.0

### 1.8.0 (2024-03-03)
- Major Update: If you experience any problems with this version, try the previous version (v1.7.2) and compare the results. If the problem persists, please report it to GitHub issues.
- The scripts in `juicer/scripts/common` used in Juicer were replaced to the original C codes and the computational time of Juicer is dramatically improved.
    - Replaced `chimeric_blacklist.awk` with `Juicer_chimeric_blacklist`
    - Replaced `dups.awk` with `Juicer_remove_duplicate`
    - Replaced `fragment.pl` with `Juicer_fragment`
    - Replaced `statistics.pl` with `Juicer_statistics`
- During this improvement I found a bug in `chimeric_blacklist.awk` where the names of non-autosomes (e.g. chrX and chrY) were not distinguished. This bug has been fixed in the new `Juicer_chimeric_blacklist`.
    - This bug seems to be related to the new version of Juicer, which uses the chromosome names chr1,chr2,chr3, not 1,2,3. I have confirmed that the output of `Juicer_chimeric_blacklist` is identical to the older Juicer.
- Added the Jupyter notebooks for the CustardPy API in the `tutorial/` directory.
- Added the notebooks to the Manual.
- Updated HiC1Dmetrics from v0.2.5 to v0.2.9

### 1.7.2 (2024-02-03)
- Modified `gethg38genome.sh` and `getmm39genome.sh` to also create the genometable file.
- Bug fix in samtools installation
- Updated SAMtools from 1.17 to 1.19.2
- Updated SRAtoolkit from 3.0.2 to 3.0.10
- Change WORKDIR from /opt to /home/ubuntu
- Installed `sudo`

### 1.7.1 (2024-01-23)
- Modified `gethg38genome.sh` and `getmm39genome.sh` to also download the refFlat data.

### 1.7.0 (2024-01-07)
- Fixed the bug in `plotHiCfeature` where TADs and loops were not being rendered
- Removed the restriction site directories of Juicer and Cooler from this repository to reduce the size
- Added `generate_binlist_from_gtfile.py` for 3DChromatin_ReplicateQC analysis
- Added `07.QualityCheck.sh` in `tutorial/Hi-C`

### 1.6.0 (2024-01-06)
- Added [Mustache](https://github.com/ay-lab/mustache)
- Added [Chromosight](https://github.com/koszullab/chromosight)
- Update Hi-C Tutorial that uses [Nakato et al.,s RPE cells](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196034) for the example data

### 1.5.0 (2023-12-18)
- AluI restriction enzyme for BAT Hi-C added to Juicer and Cooler

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
