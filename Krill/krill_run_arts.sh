#!/bin/bash
set -euo pipefail

run_conda_func () {

    CONDA_MINICONDA_PATH=~/miniconda3/etc/profile.d/conda.sh
    CONDA_ANACONDA_PATH=~/anaconda3/etc/profile.d/conda.sh

    if [ -e "$CONDA_MINICONDA_PATH" ]; then
        source "$CONDA_MINICONDA_PATH" && \
        conda activate ARTS && \
        python ~/miniconda3/envs/ARTS/arts/artspipeline1.py \
        -org METAGENOME \
        -khmms $2 \
        -t E3 \
        -rd $3 \
        -cpu $4 \
        $1 \
        ~/miniconda3/envs/ARTS/arts/reference/metagenome/
        conda deactivate
    elif [ -e "$CONDA_ANACONDA_PATH" ]; then
        source "$CONDA_ANACONDA_PATH" && \
        conda activate ARTS && \
        python ~/anaconda3/envs/ARTS/arts/artspipeline1.py \
        -org METAGENOME \
        -khmms $2 \
        -t E3 \
        -rd $3 \
        -cpu $4 \
        $1 \
        ~/anaconda3/envs/ARTS/arts/reference/metagenome/
        conda deactivate

    else
        echo "conda.sh not found in either Miniconda or Anaconda installation directories."
    fi
}

run_conda_func $1 $2 $3 $4
