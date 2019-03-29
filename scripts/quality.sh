#!/usr/bin/env bash

# Load variables
. scripts/setup.sh

mkdir -p quality

# Quality Sample data
for SAMPLE in ${SAMPLES[@]};
do
    R=reads/${SAMPLE}.fq
    mkdir -p quality/${SAMPLE}
    fastqc -o quality/${SAMPLE} $R
done
