#!/usr/bin/env bash

# Load variables
. scripts/setup.sh

mkdir -p quality

# Quality Raw Data
for SAMPLE in ${SAMPLES[@]};
do
    R=reads/${SAMPLE}.fq
    mkdir -p quality/${SAMPLE}
    fastqc -o quality/${SAMPLE} $R
done

# Quality in Alignment
#picard CollectAlignmentSummaryMetrics -Xmx2G R=refs/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa I=bam/Normal1_mappable.bam O=output.txt
