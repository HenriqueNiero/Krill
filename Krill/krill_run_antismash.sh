#!/bin/bash
set -euo pipefail

run_docker_func () {
    bash run_antismash "$1" AntiSMASH \
    --output-dir "$3" \
    --cc-mibig --rre --tigrfam \
    --skip-zip-file \
    --cb-general \
    --cb-knownclusters \
    --cb-subclusters \
    --asf \
    --pfam2go \
    --genefinding-tool prodigal-m \
    --cpus $2 
}

run_conda_func () {
    CONDA_MINICONDA_PATH=~/miniconda3/etc/profile.d/conda.sh
    CONDA_ANACONDA_PATH=~/anaconda3/etc/profile.d/conda.sh

    if [ -e "$CONDA_MINICONDA_PATH" ]; then
        source "$CONDA_MINICONDA_PATH" && \
        antismash --cc-mibig --rre --tigrfam \
        --skip-zip-file \
        --cb-general \
        --cb-knownclusters \
        --cb-subclusters \
        --asf \
        --pfam2go \
        --genefinding-tool prodigal-m \
        --cpus $2 \
        --output-dir $3 \
        $1
    elif [ -e "$CONDA_ANACONDA_PATH" ]; then
        source "$CONDA_ANACONDA_PATH" && \
        antismash --cc-mibig --rre --tigrfam \
        --skip-zip-file \
        --cb-general \
        --cb-knownclusters \
        --cb-subclusters \
        --asf \
        --pfam2go \
        --genefinding-tool prodigal-m \
        --cpus $2 \
        --output-dir $3 \
        $1
    else
        echo "conda.sh not found in either Miniconda or Anaconda installation directories."
    fi
}

if ! command -v docker &> /dev/null; then
    run_conda_func "$1" "$2" "$4"
    exit
else
    if [[ $(docker images -q antismash/standalone) ]]; then
        run_docker_func "$1" "$2" "$4"
    else
        run_conda_func "$1" "$2" "$4"
    fi
fi

