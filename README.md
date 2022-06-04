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

See [CustardPy Manual](https://custardpy.readthedocs.io/en/latest/).

## Contact

Ryuichiro Nakato: rnakato AT iqb.u-tokyo.ac.jp
