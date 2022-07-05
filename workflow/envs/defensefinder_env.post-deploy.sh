#!/usr/bin/env bash

set -euo pipefail 

# install defence finder
pip install mdmparis-defense-finder

# get database
defense-finder update
