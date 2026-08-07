#!/bin/bash
set -e

ENV=$1
INPUT=$2
OUTPUT=$3
THREADS=$4
CUTOFF=$5
PFAM=$6

echo "==================================="
echo "Running BiG-SCAPE2"
echo "Input    : $INPUT"
echo "Output   : $OUTPUT"
echo "Threads  : $THREADS"
echo "Cutoff   : $CUTOFF"
echo "==================================="

conda run -n "${ENV}" bigscape cluster \
        --input-dir "${INPUT}" \
        --output-dir "${OUTPUT}" \
        --input-mode recursive \
        --mix \
        --include-singletons \
        --gcf-cutoffs "${CUTOFF}" \
        --mibig-version 4.0 \
        --classify category \
        --pfam-path "${PFAM}" \
        --cores "${THREADS}"