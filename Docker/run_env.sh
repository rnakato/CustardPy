#!/bin/bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate $1
shift
exec "$@"
