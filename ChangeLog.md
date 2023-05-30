# Changelog

## Docker_CustardPy

### 1.1.0 (2023-05-24)
- Add `DEG_boundary_analysis` command for TAD-proximity analysis
- Add `custardpy_clustering_boundary` command for clustering TAD boundaries
- Add `custardpy_differential_DRF` command for differential DRF analysis
- Bug fix in `makeEigen.sh`
- Change name of `custardpy_mappingHiC` to `custardpy_cooler_HiC`
- Change name of `custardpy_mappingMicroC` to `custardpy_cooler_MicroC`

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

- Add tools
    - [GENOVA](https://github.com/robinweide/GENOVA)
    - [CHESS](https://chess-hic.readthedocs.io/en/latest/index.html)
    - [FIREcaller](https://github.com/yycunc/FIREcaller)
    - [STRIPENN](https://github.com/VahediLab/stripenn)
    - [HiC-Pro](https://github.com/nservant/HiC-Pro)
    - [coolpup.py](https://github.com/open2c/coolpuppy)
    - [FitHiChIP](https://ay-lab.github.io/FitHiChIP/html/index.html)
    - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

- Merged `drawTriangleRatioMulti` to `plotHiCfeature`
- Fix the error "QStandardPaths: XDG_RUNTIME_DIR points to non-existing path '/run/user/1000', please create it with 0700 permissions." when calculating PC1.
- Remove /root/.cpanm/work directory

### 0.4.1 (2023-02-24)
	- Add `libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-render-util0` libraries to avoid an error 'qt.qpa.plugin: Could not load the Qt platform plugin "xcb" in "" even though it was found.'

### 0.4.0 (2023-01-26)
	- Update R from 3.6 to 4.2
	- Add [FitHiC](https://github.com/ay-lab/fithic)
	- Add [ChIAPop](https://github.com/wh90999/ChIAPoP) for ChIA-PET analysis

### 0.3.0 (2022-10-28)
- downgrade pairtools from v1.0.1 to v0.3.0

<!--
	- change custardpy_mappingMicroC not to output BAM file that takes long time
- add `--backend cython` option to *pairtools dedup* for the consistensy to pairtools v0.3.0
-->

### 0.2.0 (2022-08-31)
- Public release
- Bug Fix for visualization
- Update Manual


### 0.1.0
- First commit

## Docker_CustardPy_Juicer (deprecated)

### 0.4.0 (2023-04-25)
- Change OS from Ubuntu16.06 to Ubuntu20.04
- Change CUDA driver from v8 to 11.0.3

### 0.3.2 (2023-02-07)
- Change the upper limit of memory in Java from `-Xmx32768m` to `-Xmx128768m`

### 0.3.1 (2022-12-12)
- Modify mega.sh to adjust CustardPy_juicer
- Add juicer_genhic.sh that generates a .hic file from the mapped data
- Deprecate `juicer_postprocessing.sh`

### 0.3.0 (2022-10-27)
- Downgrade `juicer_tools` from juicer_tools.2.13.07.jar to juicer_tools.1.22.01.jar to keep the consistency between .hic and .cool (see https://github.com/deeptools/HiCExplorer/issues/798)

### 0.2.0 (2022-08-31)
- Add restriction sites for Arima and Sau3AI
- Add restriction sites for various species
- change R installation from conda to apt

### 0.1.0
- First commit
