#!/bin/env bash
# Install latest Miniforge build for x86-64 architecture.

# Miniforge3 should take care of everything else.
curl -o ~/Miniforge3-Linux-x86_64.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash ~/Miniforge3-Linux-x86_64.sh