#!/usr/bin/env bash

# Set Alias for docker commands
alias drun="docker run -it --rm -v $(pwd):/directory -w /directory carlgira/rna-analysis:latest /bin/bash -c"

# Reference Genome
export REF=refs/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa
export GTF=refs/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.1.gtf

export SAMPLES=(SugarFed1_mappable SugarFed2_mappable Normal1_mappable Normal2_mappable Inactivated1_mappable Inactivated2_mappable Control1_mappable Control2_mappable Control3_mappable C5a1_mappable C5a2_mappable C5a3_mappable C3a1_mappable C3a2_mappable C3a3_mappable)

export LOG=log.txt
