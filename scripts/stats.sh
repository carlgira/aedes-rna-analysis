#!/usr/bin/env bash

# Load variables
. scripts/setup.sh

mkdir -p stats

DES=(c3a_vs_control c5a_vs_control inactivated_vs_normal inactivated_vs_sugarFed normal_vs_sugarFed)

for DE in ${DES[@]};
do
    echo "*** Differential Expression $DE"
    python3 src/py/filter_count.py $DE
    cat count/${DE}_norm_matrix_filter.txt | Rscript src/r/draw-heatmap.r > stats/${DE}_heatmap.pdf
done
