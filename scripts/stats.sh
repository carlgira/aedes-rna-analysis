#!/usr/bin/env bash

# Load variables
. scripts/setup.sh

# Exit this script on any error.
set -euo pipefail

mkdir -p stats

DES=(c3a_vs_control c5a_vs_control inactivated_vs_normal inactivated_vs_sugarFed normal_vs_sugarFed)

for DE in ${DES[@]};
do
    echo "Differential Expression $DE"
    drun "python3 scripts/filter_count.py $DE"
    cat count/${DE}_norm_matrix_filter.txt | Rscript scripts/draw-heatmap.r > stats/${DE}_heatmap.pdf
done
