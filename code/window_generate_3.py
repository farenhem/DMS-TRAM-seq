#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
import math
import datetime
from optparse import OptionParser
from scipy.optimize import curve_fit
from scipy.stats import rankdata

parser = OptionParser()
parser.add_option("--size", dest="SIZE", help="Integer-value size of windows to be calculated.")
parser.add_option("--path", dest="PATH", help="Full path to working directory, which contains the bedGraph and set_windows folders.")
parser.add_option("--comparison", dest="COMP", help="The prefixes of both sample conditions, connected by an underscore. Ex: condition1_condition2.")

print('Counting significant windows.')

    
(options, args) = parser.parse_args()

wkdir = str(options.PATH)
windowSize = int(options.SIZE)
comp = str(options.COMP)

os.chdir(wkdir + '/set_windows')
STEP = pd.read_csv(comp + '_' + str(windowSize) + 'bp_STEP_deltas.csv', header = 0, names = ['chr', 'start', 'end', 'ID', 'gene', 'type', 'exon', 'strand', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval', 'HepG2', 'K562', 'SG_FC', 'SG_localization'])
SLIDE = pd.read_csv(comp + '_' + str(windowSize) + 'bp_SLIDE_deltas.csv', header = 0, names = ['chr', 'start', 'end', 'ID', 'gene', 'type', 'exon', 'strand', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval', 'HepG2', 'K562', 'SG_FC', 'SG_localization'])

os.chdir(wkdir+'/reference_annot')
ID_length_name = pd.read_csv('ID_length_name.txt', sep='\t', header = 0, names = ['ID', 'gene', 'length'])
ID_type = (pd.read_csv('ensemblID_transcripttype.txt', sep='\t', names=['ID', 'type'])).fillna('')
                                                                                              
def single_annot_only(df):
    df_noNaN = df[df['ID'].notna()].reset_index(drop=True)
    df_new = df_noNaN
    for i in range(0, len(df_new)):
        if len(df_new['ID'][i].split(';')) == 1:
            continue 
        else:
            df_new.drop(index = i, inplace = True)
    df_new.reset_index(inplace=True, drop=True)
    return df_new
    
def add_annot(df):
    df_allcounts = df.merge(total_counts, on = 'ID', how='right')
    df_allcounts['sig_n'] = df_allcounts['sig_n'].fillna(value=0) 
    df_allcounts['ratio'] = df_allcounts['sig_n'] / df_allcounts['total_n']
    df_geneinfo = df_allcounts.merge(ID_length_name, on='ID', how='left')
    df_type = df_geneinfo.merge(ID_type, on='ID', how='left')
    return df_type
    
def BH_pval_correction(pvals_in):
    pvals = np.array(pvals_in)
    ranks = rankdata(pvals, method='ordinal')
    padj = np.empty(len(pvals))
    p_m_over_i = pvals * len(pvals) / ranks
    for i in range(0, len(ranks)):
        current_rank = ranks[i]
        padj[i] = min(1, min(p_m_over_i[ranks >= current_rank])) #ensures no p-val of a higher rank is lower; max of 1.
    return padj 
 
#pvalue correction over all chromosomes; rewrite csv files
STEP.dropna(subset=['p'], axis=0, inplace=True) #in RARE cases no p-val calculated (I've only seen happen once), will cause critical error in the correction below
SLIDE.dropna(subset=['p'], axis=0, inplace=True)
STEP['p_adj'] = BH_pval_correction(STEP.p.values)
SLIDE['p_adj'] = BH_pval_correction(SLIDE.p.values)

os.chdir(wkdir + '/set_windows')
STEP.to_csv(comp + '_' + str(windowSize) + 'bp_STEP_deltas_padj.csv', index = False, columns = ['chr', 'start', 'end', 'ID', 'gene', 'type', 'exon', 'strand', 'delta', 'p', 'p_adj', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval', 'HepG2', 'K562', 'SG_FC', 'SG_localization'])
SLIDE.to_csv(comp + '_' + str(windowSize) + 'bp_SLIDE_deltas_padj.csv', index = False, columns = ['chr', 'start', 'end', 'ID', 'gene', 'type', 'exon', 'strand', 'delta', 'p', 'p_adj', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval', 'HepG2', 'K562', 'SG_FC', 'SG_localization'])


single_annot_STEP = single_annot_only(STEP)
single_annot_STEP.to_csv(comp + '_' + str(windowSize) + 'bp_STEP_deltas_padj_singleannot.csv', index = False, columns = ['chr', 'start', 'end', 'ID', 'gene', 'type', 'exon', 'strand', 'delta', 'p', 'p_adj', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval', 'HepG2', 'K562', 'SG_FC', 'SG_localization'])


STEP_sig = single_annot_STEP[single_annot_STEP['p_adj'] <= 0.05]
sig_counts = pd.DataFrame(STEP_sig['ID'].value_counts(dropna=False)).reset_index()
sig_counts.rename(columns = {'ID' : 'sig_n', 'index' : 'ID'}, inplace=True)
total_counts = pd.DataFrame(single_annot_STEP['ID'].value_counts(dropna=False)).reset_index()
total_counts.rename(columns = {'ID' : 'total_n', 'index' : 'ID'}, inplace=True)

out_sig = add_annot(sig_counts)

os.chdir(wkdir+'/set_windows')
out_sig.to_csv(comp + '_' + str(windowSize) + 'bp_STEP_count_pval.csv', index = False)

print('FINISHED!')