#!/usr/bin/env bash

# Load variables
. scripts/setup.sh

mkdir -p quality

SAMPLE_FILE=$1

FASTA_FILES=($(awk 'NR>1{print $1}' $SAMPLE_FILE))

echo "*** [$(date)] [quality.sh] Fastqc $SAMPLE_FILE"
for index in ${!FASTA_FILES[*]};
 do
   FASTA=${FASTA_FILES[$index]}
   if [ ! -f $FASTA ]; then
     mkdir -p /work/quality/${FASTA}
     fastqc -o /work/quality/${FASTA} $FASTA
   fi
done

echo "*** [$(date)] [quality.sh] done"
