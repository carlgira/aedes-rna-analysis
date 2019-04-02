#!/usr/bin/env bash

# Load variables
. scripts/setup.sh

mkdir -p count

EXP_NAME=$1
SAMPLE_FILE=$2

BAM_SET0=`awk  'BEGIN { ORS=" " }; $4 == 0 {print $2}' $SAMPLE_FILE`
SIZE_SET0=`awk  '$4 == 0 {print $2}' $SAMPLE_FILE | wc -l`

BAM_SET1=`awk  'BEGIN { ORS=" " }; $4 == 1 {print $2}' $SAMPLE_FILE`
SIZE_SET1=`awk '$4 == 1 {print $2}' $SAMPLE_FILE | wc -l`

COUNT_FILE=count/${EXP_NAME}_count.tsv
MATRIX_FILE=count/${EXP_NAME}_matrix.tsv
EDGER_FILE=count/${EXP_NAME}_edger.tsv

# Generate the counts.
echo "*** Counting features with: $GTF"
featureCounts -a $GTF -g gene_id -o $COUNT_FILE $BAM_SET0 $BAM_SET1

# Remove header and other columns
#tail -n +2 $COUNT_FILE | cut -f1,7- > temp.txt
#mv temp.txt $COUNT_FILE

# Run the EdgeR method
echo "*** Running EdgeR."
cut -f1,7- $COUNT_FILE | Rscript src/r/edger.r ${SIZE_SET0}x${SIZE_SET1} $MATRIX_FILE > $EDGER_FILE
