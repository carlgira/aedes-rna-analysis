
. scripts/setup.sh

EXP_NAMES=(c3a_vs_control c5a_vs_control normal_vs_sugarfed inactivated_vs_sugarfed inactivated_vs_normal)
SAMPLE_FILES=(setup/analysis_samples_vitro_exp1.tsv setup/analysis_samples_vitro_exp2.tsv setup/analysis_samples_vivo_exp1.tsv setup/analysis_samples_vivo_exp2.tsv setup/analysis_samples_vivo_exp3.tsv)
ALIGNER='star'
NCPU=2

# Prepare data
echo "*** [$(date)] [run-all.sh] Prepare data"
. scripts/prepare-data.sh

# Align and count all
echo "*** [$(date)] [run-all.sh] Aligning all samples"
for index in ${!EXP_NAMES[*]};
 do
   EXP_NAME=${EXP_NAMES[$index]}
   SAMPLE_FILE=${SAMPLE_FILES[$index]}
  . scripts/align-and-count.sh $EXP_NAME $SAMPLE_FILE
done

echo "*** [$(date)] [run-all.sh] Done"
