#!/usr/bin/env bash

set -euo pipefail 

# install hammer
mamba install -c bioconda hmmer -y

# install defence finder
pip3 install mdmparis-defense-finder==1.0.8 biopython

# get database
# defense-finder update
macsydata install -U -u --org mdmparis defense-finder-models
