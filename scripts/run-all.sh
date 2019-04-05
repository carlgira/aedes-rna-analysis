
. scripts/setup.sh

EXP_NAMES=(c3a_vs_control c5a_vs_control normal_vs_sugarfed inactivated_vs_sugarfed inactivated_vs_normal)
ALIGNER='star'
NCPU=8
CPM=0.5
PVALUE=0.05

# Prepare data
echo "*** [$(date)] [run-all.sh] Prepare data"
#. scripts/prepare-data.sh

# Align and count all
echo "*** [$(date)] [run-all.sh] Aligning all samples"
for index in ${!EXP_NAMES[*]};
 do
   EXP_NAME=${EXP_NAMES[$index]}
   SAMPLE_FILE=/work/setup/${EXP_NAME}.tsv
  #. scripts/align-and-count.sh $EXP_NAME $SAMPLE_FILE
done

# Fastqc fasta files
echo "*** [$(date)] [run-all.sh] Fastqc samples"
for index in ${!EXP_NAMES[*]};
 do
   SAMPLE_FILE=${SAMPLE_FILES[$index]}
  #. scripts/quality.sh $SAMPLE_FILE
done

# MultiQC Report
#multiqc -o output
# Create R reports
echo "*** [$(date)] [run-all.sh] R Report"
for index in ${!EXP_NAMES[*]};
 do
   EXP_NAME=${EXP_NAMES[$index]}
   R -e "rmarkdown::render('de-template.Rmd',output_file='output/${EXP_NAME}.html', params=list(expName='${EXP_NAME}', cpm=${CPM}, pvalue=${PVALUE}))"
done

echo "*** [$(date)] [run-all.sh] Done"
