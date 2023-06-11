#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate ReplicateQC
cd /opt/3DChromatin_ReplicateQC

/opt/conda/envs/ReplicateQC/bin/3DChromatin_ReplicateQC $@
