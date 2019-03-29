#!/usr/bin/env bash

# Load variables
. scripts/setup.sh

# Create output folder
mkdir -p bam

# Iterate over each sample
for SAMPLE in ${SAMPLES[@]};
do
    R=reads/${SAMPLE}.fq
    BAM=bam/${SAMPLE}.bam

    if [ ! -f $BAM ]; then
      echo "*** Aligning: $BAM $(date)"
      hisat2 $REF -U $R | samtools sort > $BAM
      samtools index $BAM
    fi
done
