# CustardPy: Docker image for 3D genome analysis

<img src = "image/CustardPy.jpg" width = 700ptx>

## Major Release! (version 1)

This repository contains
- Source code of **CustardPy** (PyPI),
- Dockerfile of **CustardPy** Docker image, 
- [Full Manual](https://custardpy.readthedocs.io), and
- Tutorial data of Hi-C and Micro-C analysis using demo data.

## 0. Changelog

See [Changelog](https://github.com/rnakato/CustardPy/blob/main/ChangeLog.md)

## 1. Installation

Docker image is available at [DockerHub](https://hub.docker.com/r/rnakato/custardpy).

### 1.1 Docker

To use docker command, type:

    # pull docker image
    docker pull rnakato/custardpy

    # container login
    docker run --rm -it rnakato/custardpy /bin/bash
    # execute a command
    docker run --rm -it -v (your directory):/opt/work rnakato/custardpy <command>

When calling loops using Juicer HICCUPS, supply the ``--gpus all`` option to allow GPU computation (GPU card needed):

    docker run --gpus all -it --rm -it -v (your directory):/opt/work rnakato/custardpy call_HiCCUPS.sh

- user:password
    - ubuntu:ubuntu

### 1.2 Singularity

Singularity can also be used to execute the docker image.
You can build the singularity file (.sif) of CustardPy with this command:

    singularity build custardpy.sif docker://rnakato/custardpy
    
Instead, you can download the CustardPy singularity image from our [Dropbox](https://www.dropbox.com/scl/fo/lptb68dirr9wcncy77wsv/h?rlkey=whhcaxuvxd1cz4fqoeyzy63bf&dl=0).

Then you can run CustardPy with the command:

    singularity exec custardpy.sif <command>

Singularity will automatically mount the current directory. If you want to access the files in the other directory, use the `--bind` option, for instance:

    singularity exec --bind /work custardpy.sif <command>

This command mounts `/work` directory.

When calling loops using Juicer HICCUPS, supply ``--nv`` option to allow GPU computation (GPU card needed):

    singularity exec --bind /work custardpy.sif call_HiCCUPS.sh

## 2. Quickstart

    # download Churros/tutorial directory
    git clone https://github.com/rnakato/CustardPy.git
    cd CustardPy/tutorial/Hi-C/

    # download fastq and genome data and make index
    bash 00_getdata.sh

    # Execute Juicer pipeline
    bash QuickStart_juicer.sh

## 3. Usage

See https://custardpy.readthedocs.io for the detailed Manual.

## 4. Build Docker image from Dockerfile

First clone and move to the repository

    git clone https://github.com/rnakato/CustardPy.git
    cd CustardPy/Docker

Then type:

    docker build -f Dokerfile.<version> -t <account>/custardpy_juicer .

## 6. Reference

- Nakato R, Sakata T, Wang J, Nagai LAE, Nagaoka Y, Oba GM, Bando M, Shirahige K, Context-dependent perturbations in chromatin folding and the transcriptome by cohesin and related factors, *Nature Communications*, 2023. doi: [10.1038/s41467-023-41316-4](https://www.nature.com/articles/s41467-023-41316-4)
