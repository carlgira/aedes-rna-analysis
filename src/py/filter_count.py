import pandas as pd
import sys
file = sys.argv[1]
df = pd.read_csv('count/' + file + '_results_edger.txt', delimiter='\t')
norm = pd.read_csv('count/' + file + '_norm_matrix.txt', delimiter='\t', index_col='id')
norm[norm.index.isin(df[(df.PValue < 0.01) & (df.FDR < 0.05)].index)].to_csv('count/' + file + '_norm_matrix_filter.txt', sep='\t')
