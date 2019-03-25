#!/usr/bin/env bash

# Load variables
. scripts/setup.sh

set -euxo pipefail

mkdir -p quality

# Quality Sample data
for SAMPLE in ${SAMPLES[@]};
do
    R=reads/${SAMPLE}.fq
    mkdir -p quality/${SAMPLE}
    fastqc -o quality/${SAMPLE} $R
done
