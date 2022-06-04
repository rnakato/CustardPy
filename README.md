# CustardPy

Hi-C analysis tools by Python3 and Docker


## Requirements

The following are required before installing CustardPy:

- Python 3.7+

## Installation

### From PyPI

Core component of CustardPy can by installed using pip:

    pip3 install custardpy

### Docker image

Docker image of CustardPy that contains additional scripts for Hi-C/Micro-C analysis is available at [DockerHub](https://hub.docker.com/repository/docker/rnakato/custardpy/general).

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

See [Custardpy Manual](https://custardpy-juicer.readthedocs.io/en/latest/index.html).
