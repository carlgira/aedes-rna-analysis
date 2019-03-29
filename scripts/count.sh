#!/usr/bin/env bash

# Load variables
. scripts/setup.sh

mkdir -p count

# Generate the counts.
echo "*** Counting features with: $GTF"
featureCounts -a $GTF -g gene_id -o count/vivo_counts.txt  bam/Normal*.bam bam/SugarFed*.bam bam/Inactivated*.bam
featureCounts -a $GTF -g gene_id -o count/vitro_counts.txt  bam/Control*.bam bam/C5a*.bam bam/C3a*.bam

# Separate InVivo Counts
cat count/vivo_counts.txt | cut -f 1,7-10 > count/normal_vs_sugarFed_counts.txt
cat count/vivo_counts.txt | awk '{print($1 "\t" $11 "\t" $12 "\t" $9 "\t" $10)}' > count/inactivated_vs_sugarFed_counts.txt
cat count/vivo_counts.txt | awk '{print($1 "\t" $11 "\t" $12 "\t" $7 "\t" $8)}' > count/inactivated_vs_normal_counts.txt

# Separate InVitro Counts
cat count/vitro_counts.txt | awk '{print($1 "\t" $13 "\t" $14 "\t" $15 "\t" $7 "\t" $8 "\t" $9)}' > count/c3a_vs_control_counts.txt
cat count/vitro_counts.txt | awk '{print($1 "\t" $10 "\t" $11 "\t" $12 "\t" $7 "\t" $8 "\t" $9)}' > count/c5a_vs_control_counts.txt

# Run the EdgeR method
echo "*** Running EdgeR."
cat count/normal_vs_sugarFed_counts.txt | Rscript src/r/edger.r 2x2 count/normal_vs_sugarFed_norm_matrix.txt > count/normal_vs_sugarFed_results_edger.txt
cat count/inactivated_vs_sugarFed_counts.txt | Rscript src/r/edger.r 2x2 count/inactivated_vs_sugarFed_norm_matrix.txt > count/inactivated_vs_sugarFed_results_edger.txt
cat count/inactivated_vs_normal_counts.txt | Rscript src/r/edger.r 2x2 count/inactivated_vs_normal_norm_matrix.txt > count/inactivated_vs_normal_results_edger.txt

cat count/c3a_vs_control_counts.txt | Rscript src/r/edger.r 3x3 count/c3a_vs_control_norm_matrix.txt > count/c3a_vs_control_results_edger.txt
cat count/c5a_vs_control_counts.txt | Rscript src/r/edger.r 3x3 count/c5a_vs_control_norm_matrix.txt > count/c5a_vs_control_results_edger.txt
