# CustardPy

Hi-C analysis tools by Python3 and Docker


## Requirements

The following are required before installing CustardPy:

- Python 3.7+

## Installation

### From PyPI

Core components of CustardPy can by installed using pip:

    pip3 install custardpy

### Docker image

We recommend to use the [CustardPy Docker image](https://hub.docker.com/r/rnakato/custardpy) that contains additional scripts for Hi-C/Micro-C analysis.

#### Docker 
To use docker command, type:

    docker pull rnakato/custardpy
    docker run -it --rm rnakato/custardpy <command>

#### Singularity

Singularity can also be used to execute the docker image:

    singularity build custardpy.sif docker://rnakato/custardpy
    singularity exec custardpy.sif <command>

Singularity mounts the current directory automatically. If you access the files in the other directory, please mount by `--bind` option, for instance:

    singularity exec --bind /work custardpy.sif <command>
    
This command mounts `/work` directory.

## Usage

See https://custardpy.readthedocs.io for the detailed Manual.

## Singularity images

Prebuild singularity images (version 3.8.5) are available on our [Google Drive](https://drive.google.com/drive/folders/1wkw19qPKm8lnWorBXu937zOvwIANbiGw?usp=sharing).

## Reference

- Nakato R, Sakata T, Wang J, Nagai LAE, Oba GM, Bando M, Shirahige K, Context-dependent 3D genome regulation by cohesin and related factors, bioRxiv, 2022. doi: [10.1101/2022.05.24.493188](https://www.biorxiv.org/content/10.1101/2022.05.24.493188v1)
