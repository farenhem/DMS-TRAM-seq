import os
import pandas as pd
import numpy as np
from scipy.stats import rankdata
import argparse

def BH_pval_correction(pvals_in):
        pvals = np.array(pvals_in)
        ranks = rankdata(pvals, method='ordinal')
        padj = np.empty(len(pvals))
        p_m_over_i = pvals * len(pvals) / ranks
        for i in range(0, len(ranks)):
            current_rank = ranks[i]
            padj[i] = min(1, min(p_m_over_i[ranks >= current_rank])) #ensures no p-val of a higher rank is lower; max of 1
        return padj
        
parser = argparse.ArgumentParser(description="Specify the input file.")
parser.add_argument('infile', type=str, help='Path to the input file.')
args = parser.parse_args()   

infile = str(args.infile)

filetype = infile.split('.')[-1] 
if filetype == 'txt' or filetype == 'bed':
    sep = '\t'
elif filetype == 'csv':
    sep = ','

windows = pd.read_csv(infile, index_col=False, sep=sep)
windows['padj'] = BH_pval_correction(windows['p'])

outfile = '.'.join(infile.split('.')[:-1]) + '_padj'

windows.to_csv(outfile + '.csv', index=False, columns = ['gene', 'ID', 'type', 'chr', 'strand', 'region', 'exon', 'start', 'end',
                                                         'x', 'delta', 'p', 'padj', 'r_con1', 'r_con2', 'r_between', 'con1_avg', 'con2_avg', 'avg_pval', 
                                                         'con1_gini', 'con2_gini', 'gini_pval', 'n', 'len'])
windows.to_csv(outfile + '.txt', index=False, sep='\t', columns = ['gene', 'ID', 'type', 'chr', 'strand', 'region', 'exon', 'start', 'end',
                                                                    'x', 'delta', 'p', 'padj', 'r_con1', 'r_con2', 'r_between', 'con1_avg', 'con2_avg', 'avg_pval', 
                                                                    'con1_gini', 'con2_gini', 'gini_pval', 'n', 'len'])