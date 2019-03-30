#!/usr/bin/env bash


. scripts/setup.sh

# Create output folder
mkdir -p refs

# VectorBase - Aedes Reference Genome
if [ ! -f $REF ]; then
  echo "*** Get Genome $(date)"
  curl -L https://www.vectorbase.org/download/aedes-aegypti-lvpagwgchromosomesaaegl5fagz -o refs/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa.gz
  gunzip ${REF}.gz

  echo "*** Index Genome $(date)"
  hisat2-build $REF $REF
fi

# VectorBase - Aedes Gene Feature File
if [ ! -f $GTF ]; then
  echo "*** Get GTF $(date)"
  curl -L https://www.vectorbase.org/download/aedes-aegypti-lvpagwgbasefeaturesaaegl51gtfgz -o refs/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.1.gtf.gz
  gunzip ${GTF}.gz
fi

# Unzip Sample data
echo "*** Unzip sample data"
for SAMPLE in ${SAMPLES[@]};
do
    R=reads/${SAMPLE}.fq
    GZ=reads/${SAMPLE}.fq.gz

    if [ ! -f $R ]; then
      gunzip $GZ
    fi
done
