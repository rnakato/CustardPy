#!/bin/bash
export MAMBA_ROOT_PREFIX=/opt/micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate $1
shift
exec "$@"
