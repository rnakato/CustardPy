FROM rnakato/ubuntu_gpu_20.04:cuda11.0.3-cudnn8 AS base
LABEL maintainer="Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>"

# For sorting, LC_ALL is C
ENV LC_ALL=C
ENV NVIDIA_VISIBLE_DEVICES=all
ENV NVIDIA_DRIVER_CAPABILITIES=all

WORKDIR /opt/
USER root

SHELL ["/bin/bash", "-c"]

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    apt-utils \
    automake \
    bc \
    build-essential \
    bzip2 \
    default-jdk \
    clang \
    cmake \
    curl \
    fastqc \
    ffmpeg \
    gawk \
    gcc \
    git \
    gpg-agent \
    imagemagick \
    less \
    libbz2-dev \
    libclang-dev \
    libcurl4-gnutls-dev \
    libfontconfig1-dev \
    libfreetype-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgit2-dev \
    libglpk-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    liblz4-tool \
    libncurses-dev \
    libncurses5 \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libuv1-dev \
    libwebp-dev \
    libxkbcommon-x11-0 \
    libxcb-icccm4 \
    libxcb-image0 \
    libxcb-keysyms1 \
    libxcb-render-util0 \
    libxkbcommon-x11-0 \
    libxml2-dev \
    libz-dev \
    locales \
    make \
    pigz \
    qtcreator \
    unzip \
    zlib1g-dev \
    && echo "deb https://cran.rstudio.com/bin/linux/ubuntu focal-cran40/" | tee -a /etc/apt/sources.list \
    && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 51716619E084DAB9 \
    && apt-get update \
    && apt-get install -y --no-install-recommends r-base r-base-core r-recommended r-base-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# R packages
ENV Ncpus=16
RUN R -e "install.packages('BiocManager')" \
    && R -e "BiocManager::install(version = '3.22', ask = FALSE)" \
    && R CMD javareconf \
    && R -e "install.packages(c('data.table', \
                                'devtools', \
                                'dplyr', \
                                'fdrtool', \
                                'ggplot2', \
                                'ggpubr', \
                                'ggsci', \
                                'hash', \
                                'Nozzle.R1', \
                                'optparse', \
                                'parallel', \
                                'plyr', \
                                'Rcpp', \
			                 	'remotes', \
                                'reshape2', \
                                'sm', \
                                'splines', \
                                'strawr', \
                                'stringr', \
                                'tools', \
                                'tidyr', \
                                'tidyverse'))" \
    && R -e "remotes::install_github('hadley/lazyeval')" \
    && R -e "BiocManager::install(c('BiocGenerics', \
                                    'Biostrings', \
                                    'edgeR', \
                                    'GenomeInfoDb', \
                                    'GenomicAlignments', \
                                    'GenomicRanges', \
                                    'IRanges', \
                                    'matrixStats', \
                                    'rhdf5', \
                                    'rtracklayer', \
                                    'S4Vectors', \
                                    'ShortRead'))"

# CALDER2: https://github.com/CSOgroup/CALDER2

RUN wget --quiet https://x.gd/88PPP -O CALDER2.tar.gz \
    && R -e " install.packages(c('R.utils', 'doParallel', 'ape', 'dendextend', 'fitdistrplus', 'igraph', 'rARPACK', 'factoextra', 'fields'))" \
    && R CMD INSTALL /opt/CALDER2.tar.gz \
    && rm -rf /opt/CALDER2.tar.gz

RUN R -e 'remotes::install_github(c("robinweide/GENOVA", "yycunc/FIREcaller"))'
# R -e 'remotes::install_url("https://github.com/SooLee/plotosaurus/archive/0.9.2.zip")' \
#    && R -e 'remotes::install_bitbucket("chicagoTeam/Chicago", subdir="Chicago")' \
#    && R -e 'remotes::install_bitbucket("chicagoTeam/Chicago", subdir="PCHiCdata")'

COPY hicrep_1.12.2.tar.gz hicrep_1.12.2.tar.gz
RUN R -e 'install.packages("/opt/hicrep_1.12.2.tar.gz", repo = NULL, type = "source")' \
    && rm hicrep_1.12.2.tar.gz

