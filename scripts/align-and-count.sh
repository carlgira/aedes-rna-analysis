# Load variables
. scripts/setup.sh

# Output folders
mkdir -p count
mkdir -p bam

if [ $# -ne 2 ]
  then
    echo "Parammeters (EXP_NAME -> Experiment name, SAMPLE_FILE -> Sample conf file)"
    exit 1
fi

# Parammeters (EXP_NAME -> Experiment name, SAMPLE_FILE -> Sample conf file)
EXP_NAME=$1
SAMPLE_FILE=$2

FASTA_FILES=($(awk 'NR>1{print $1}' $SAMPLE_FILE))
BAM_FILES=($(awk 'NR>1{print $2}' $SAMPLE_FILE))
COUNT_FILE=count/${EXP_NAME}_count.tsv

# Aligning fasta files
echo "*** [$(date)] [align-and-count.sh] Aligning fasta files $EXP_NAME"
for index in ${!FASTA_FILES[*]};
 do
   FASTA=${FASTA_FILES[$index]}
   BAM=${BAM_FILES[$index]}

   if [ ! -f $BAM ]; then
     echo "*** Aligning: $BAM $(date)"

     if [ $ALIGNER = "star" ]; then
       STAR --genomeDir /work/refs --runThreadN $NCPU \
       --readFilesIn $FASTA \
       --outSAMtype BAM SortedByCoordinate \
       --limitGenomeGenerateRAM 26000000000 \
       --outStd BAM_SortedByCoordinate > $BAM
     else
       hisat2 $REF -U $FASTA | samtools sort > $BAM
     fi
     samtools index $BAM
   fi
done

BAM_FILES=`awk 'BEGIN {ORS=" "}; NR>1{print $2}' $SAMPLE_FILE`

# Generate counts.
echo "*** [$(date)] [align-and-count.sh] Counting features"
if [ ! -f $COUNT_FILE ]; then
    featureCounts -a $GTF -g gene_id -o $COUNT_FILE $BAM_FILES
fi

echo "*** [$(date)] [align-and-count.sh] Done"
