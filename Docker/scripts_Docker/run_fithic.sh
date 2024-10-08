#!/bin/bash
export MAMBA_ROOT_PREFIX=/opt/micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate hic-pro
echo "bash /opt/FitHiChIP/$@"
bash /opt/FitHiChIP/$@
