#!/usr/bin/env bash

# Reference Genome
export REF=refs/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa

# GTF
export GTF=refs/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.1.gtf

# List of samples
export SAMPLES=(SugarFed1_mappable SugarFed2_mappable Normal1_mappable Normal2_mappable Inactivated1_mappable Inactivated2_mappable Control1_mappable Control2_mappable Control3_mappable C5a1_mappable C5a2_mappable C5a3_mappable C3a1_mappable C3a2_mappable C3a3_mappable)

# Log file
export LOG=log.txt

set -euxo pipefail
