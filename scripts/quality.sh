#!/usr/bin/env bash

# Load variables
. scripts/setup.sh

mkdir -p quality

SAMPLE_FILE=$1

FASTA_FILES=($(awk 'NR>1{print $1}' $SAMPLE_FILE))
SAMPLE_NAMES=($(awk 'NR>1{print $3}' $SAMPLE_FILE))

echo "*** [$(date)] [quality.sh] Fastqc $SAMPLE_FILE"
for index in ${!FASTA_FILES[*]};
 do
   FASTA=${FASTA_FILES[$index]}
   SAMPLE_NAME=${SAMPLE_NAMES[$index]}
   Q_DIR=/work/quality/${SAMPLE_NAME}

   if [ ! -f $Q_DIR ]; then
     mkdir -p $Q_DIR
     fastqc -o /work/quality/${SAMPLE_NAME} $FASTA
   fi
done

echo "*** [$(date)] [quality.sh] done"
