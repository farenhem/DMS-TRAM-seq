#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
from optparse import OptionParser
from scipy.stats import pearsonr
from scipy import stats
from scipy.stats import ttest_ind

parser = OptionParser()
parser.add_option("--chr", dest="CHR", help="Chromosome name. Ex: chr1")
parser.add_option("--bgpath", dest="BGPATH", help="Path to bedGraph files for all samples.")
parser.add_option("--ref", dest="REFPATH", help="Path to reference annotation file (bed) with canonical exons for all transcripts.") 
parser.add_option("--out", dest="OUTPATH", help='Destination of path for output files, also containing the temp folders "gene_annot" and "gene_bedGraph".')
parser.add_option("--winsize", dest="SIZE", help="Integer value for window size (number of valid datapoints per window.")

(options, args) = parser.parse_args()
chrom_in = str(options.CHR)
bgpath = str(options.BGPATH)
outpath = str(options.OUTPATH)
ref_file = str(options.REFPATH)
windowSize = int(options.SIZE)


def winsor_scale(col):
        arr = np.array(col)
        win_arr = stats.mstats.winsorize(arr, limits = [0.00, 0.05]).data #winsorize with limits of 0%, 95%
        norm_arr = (win_arr - win_arr.min())/(win_arr.max() - win_arr.min()) #scale to 1
        return norm_arr

def gini (rate_series):
    rate_list = list(rate_series)
    n = len(rate_list)
    denom = 2 * n * np.sum(rate_list)
    numer = 0
    if denom == 0:
        gini = 0
    else:
        for i in range(n):
            for j in range(n):
                x_i, x_j = rate_list[i], rate_list[j]
                value = abs(x_i - x_j)
                numer += value
        gini = numer / denom
    return gini
    
def window_calc(df): #input df will already be sliced to only the window of interest
    #find delta values between the two conditions, paired by replicate
    if (max(df.sample1.values) == 0) | (max(df.sample2.values) == 0) | (max(df.sample3.values) == 0) | \
    (max(df.sample4.values) == 0) | (max(df.sample5.values) == 0) | (max(df.sample6.values) == 0): #skips windows where any of the samples contain all 0's. Deleted later with dropna
        delta_mean = np.nan
        p = np.nan
        r_con1 = np.nan
        r_con2 = np.nan
        r_between = np.nan
        con1_avg = np.nan
        con2_avg = np.nan
        avg_pval = np.nan
        con1_gini = np.nan
        con2_gini = np.nan
        gini_pval = np.nan
    elif len(df) < 0.5 * windowSize: #windows that are less than half the size of the intended window length are skipped for now; will be merged with the prior window and recalculated later
        delta_mean = np.nan
        p = np.nan
        r_con1 = np.nan
        r_con2 = np.nan
        r_between = np.nan
        con1_avg = np.nan
        con2_avg = np.nan
        avg_pval = np.nan
        con1_gini = np.nan
        con2_gini = np.nan
        gini_pval = np.nan
    else:
        diffs_1 = []
        diffs_2 = []
        diffs_3 = []
        for x in range (0, len(df)):
            diffs_1.append(abs(df.sample4.values[x] - df.sample1.values[x]))
            diffs_2.append(abs(df.sample5.values[x] - df.sample2.values[x]))
            diffs_3.append(abs(df.sample6.values[x] - df.sample3.values[x]))
        delta1 = sum(diffs_1) / len(diffs_1)
        delta2 = sum(diffs_2) / len(diffs_2)
        delta3 = sum(diffs_3) / len(diffs_3)
        delta_mean = np.mean([delta1, delta2, delta3])
        
        #calculate p-values by comparing pearson's R for replicates of the same condition vs between conditions

        rwithin = []
        rbetween = []
        
        rwithin.append(pearsonr(df.sample1.values, df.sample2.values).statistic)
        rwithin.append(pearsonr(df.sample1.values, df.sample3.values).statistic)
        rwithin.append(pearsonr(df.sample2.values, df.sample3.values).statistic)
        rwithin.append(pearsonr(df.sample4.values, df.sample5.values).statistic)
        rwithin.append(pearsonr(df.sample4.values, df.sample6.values).statistic)
        rwithin.append(pearsonr(df.sample5.values, df.sample6.values).statistic)
        
        r_con1 = np.average(rwithin[0:3])
        r_con2 = np.average(rwithin[3:6])
        
        rbetween.append(pearsonr(df.sample1.values, df.sample4.values).statistic)
        rbetween.append(pearsonr(df.sample1.values, df.sample5.values).statistic)
        rbetween.append(pearsonr(df.sample1.values, df.sample6.values).statistic)
        rbetween.append(pearsonr(df.sample2.values, df.sample4.values).statistic)
        rbetween.append(pearsonr(df.sample2.values, df.sample5.values).statistic)
        rbetween.append(pearsonr(df.sample2.values, df.sample6.values).statistic)
        rbetween.append(pearsonr(df.sample3.values, df.sample4.values).statistic)
        rbetween.append(pearsonr(df.sample3.values, df.sample5.values).statistic)
        rbetween.append(pearsonr(df.sample3.values, df.sample6.values).statistic)
        
        r_between = np.average([rbetween[0], rbetween[4], rbetween[8]])

        p = ttest_ind(rbetween,rwithin, equal_var=False, alternative='less').pvalue
        
        sample1_avg = sum(df.sample1.values) / len(df.sample1.values)
        sample2_avg = sum(df.sample2.values) / len(df.sample2.values)
        sample3_avg = sum(df.sample3.values) / len(df.sample3.values)
        sample4_avg = sum(df.sample4.values) / len(df.sample4.values)
        sample5_avg = sum(df.sample5.values) / len(df.sample5.values)
        sample6_avg = sum(df.sample6.values) / len(df.sample6.values)
        
        con1_avg = (sample1_avg + sample2_avg + sample3_avg) / 3
        con2_avg = (sample4_avg + sample5_avg + sample6_avg) / 3
        avg_pval = ttest_ind([sample1_avg, sample2_avg, sample3_avg], [sample4_avg, sample5_avg, sample6_avg]).pvalue
        
        sample1_gini = gini(stats.mstats.winsorize(df.sample1.values, limits=[0, 0.05]).data)
        sample2_gini = gini(stats.mstats.winsorize(df.sample2.values, limits=[0, 0.05]).data)
        sample3_gini = gini(stats.mstats.winsorize(df.sample3.values, limits=[0, 0.05]).data)
        sample4_gini = gini(stats.mstats.winsorize(df.sample4.values, limits=[0, 0.05]).data)
        sample5_gini = gini(stats.mstats.winsorize(df.sample5.values, limits=[0, 0.05]).data)
        sample6_gini = gini(stats.mstats.winsorize(df.sample6.values, limits=[0, 0.05]).data)
        
        con1_gini = (sample1_gini + sample2_gini + sample3_gini) / 3
        con2_gini = (sample4_gini + sample5_gini + sample6_gini) / 3
        gini_pval = ttest_ind([sample1_gini, sample2_gini, sample3_gini], [sample4_gini, sample5_gini, sample6_gini]).pvalue
        
    return delta_mean, p, r_con1, r_con2, r_between, con1_avg, con2_avg, avg_pval, con1_gini, con2_gini, gini_pval
    

annot = pd.read_csv(ref_file, sep='\t', header=None, names=['chr', 'start','end', 'region', 'exon', 'name', 'strand', 'type', 'ID'])
os.chdir(outpath)

ID_list = annot[annot['chr'] == chrom_in].dropna(subset=['name']).ID.unique()

step_output = pd.DataFrame(columns=['gene', 'ID', 'type', 'chr', 'strand', 'region', 'exon', 'start', 'end',
                               'x', 'delta', 'p', 'r_con1', 'r_con2', 'r_between',
                               'con1_avg', 'con2_avg', 'avg_pval', 
                               'con1_gini', 'con2_gini', 'gini_pval', 'n', 'len'])
                               
slide_output = pd.DataFrame(columns=['gene', 'ID', 'type', 'chr', 'strand', 'region', 'exon', 'start', 'end',
                               'x', 'delta', 'p', 'r_con1', 'r_con2', 'r_between',
                               'con1_avg', 'con2_avg', 'avg_pval', 
                               'con1_gini', 'con2_gini', 'gini_pval', 'n', 'len'])

for ID_in in ID_list: 
    gene_annot = annot[annot['ID'] == ID_in]
    gene_name = gene_annot.name.unique()[0]
    gene_type = gene_annot.type.unique()[0]
    
    if gene_annot.strand.unique()[0] == '+':
        strand = '+'
        strand_name = 'pos'
        gene_annot = gene_annot.sort_values(by='start', ascending=True)
        
    elif gene_annot.strand.unique()[0] == '-':
        strand = '-'
        strand_name = 'neg'
        gene_annot = gene_annot.sort_values(by='start', ascending=False)
        
    gene_annot.reset_index(drop=True, inplace=True)

    (gene_annot.sort_values(by='start', ascending=True)).to_csv('gene_annot/' + ID_in + '_exons.bed', sep='\t', header=False, index=False,
                     columns=['chr', 'start','end', 'region', 'exon'])
                     
    if len(gene_annot.chr.unique()) > 1: #skip transcripts encoded by multiple DNA regions on multiple chromosomes
        continue
    if len(gene_annot.ID.unique()) > 1: #skip transcripts encoded by multiple DNA regions
        continue
        
    os.system('bedtools intersect -a ' + bgpath + '/allSamples_' + chrom_in + '_mmRate_' + strand_name + '.bedGraph -b gene_annot/' + ID_in + '_exons.bed > gene_bedGraph/' + ID_in  + '.bedGraph')
    os.system("bedmap --echo --echo-map-id-uniq --echo-map-score --delim '\t' gene_bedGraph/" + ID_in + '.bedGraph gene_annot/' + ID_in + '_exons.bed > gene_bedGraph/' + ID_in + '_annot.bedGraph')  
    gene = pd.read_csv('gene_bedGraph/' + ID_in + '_annot.bedGraph', sep='\t', header=None, names=['chr', 'start', 'end', 'sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6', 'region', 'exon'])

    if gene_type in ['snoRNA', 'misc_RNA', 'snRNA', 'scaRNA', 'scRNA', 'miRNA'] and len(gene) < 25:  #skip small RNAs with less than 25 valid bases
        continue
    elif gene_type not in ['snoRNA', 'misc_RNA', 'snRNA', 'scaRNA', 'scRNA', 'miRNA'] and len(gene) < 100: #all other gene types (lncRNA, mRNA, etc) need to have at least 100 valid bases 
        continue

    if strand == '+':
        gene.sort_values(by='end', ascending=True, inplace=True) #should already be sorted this way
    elif strand == '-':
        gene.sort_values(by='end', ascending=False, inplace=True) #reverses df for neg stranded genes
        
    gene.reset_index(drop=True, inplace=True)
    
    for sample in ['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6']:
        gene[sample+'_wins'] = winsor_scale(gene[sample])
    
    
    data = {}
    binSize = windowSize - 1 #STEP_WINDOWS: EACH WINDOW STARTS THE BASE AFTER THE PRIOR ONE ENDS
    for i in range(0, len(gene), windowSize):
        delta_mean, p, r_con1, r_con2, r_between, con1_avg, con2_avg, avg_pval, con1_gini, con2_gini, gini_pval = window_calc(gene.loc[i:i+binSize, :])
        data[len(data)] = {'start' : min(gene.loc[i:i+binSize, :]['end']), #in gene df, 'end' coord is index1 inclusive position
             'end' : max(gene.loc[i:i+binSize, :]['end']),
             'region': ','.join(gene.loc[i:i+binSize, :].region.astype(str).unique()),
             'exon' : ','.join(gene.loc[i:i+binSize, :].exon.astype(int).astype(str).unique()),
             'delta' : delta_mean, 'p' : p, 
             'r_con1' : r_con1, 'r_con2' : r_con2, 'r_between' : r_between,
             'con1_avg' : con1_avg, 'con2_avg' : con2_avg, 'avg_pval' : avg_pval, 
             'con1_gini' : con1_gini, 'con2_gini' : con2_gini, 'gini_pval' : gini_pval,
             'n': len(gene.loc[i:i+binSize, :])}
              
    step_windows = pd.DataFrame.from_dict(data=data, orient = 'index', columns=['gene', 'ID', 'type', 'chr', 'strand', 'region', 'exon', 'start', 'end',
                                                                           'x', 'delta', 'p', 'r_con1', 'r_con2', 'r_between', 
                                                                           'con1_avg', 'con2_avg', 'avg_pval', 
                                                                           'con1_gini', 'con2_gini', 'gini_pval', 'n', 'len'])
    step_windows['gene'] = gene_name
    step_windows['ID'] = gene_annot.ID.unique()[0]
    step_windows['type'] = gene_annot.type.unique()[0]
    step_windows['chr'] = gene_annot.chr.unique()[0] 
    step_windows['strand'] = gene_annot.strand.unique()[0]
    step_windows['err'] = ''

    step_windows['x'] = range(1, len(step_windows) + 1)

    if step_windows.loc[len(step_windows)-1].n < (0.5 * windowSize):
        #when the last window is too small (less than half the normal size), combine the last two windows & recalculate
        new_start = min(step_windows.tail(2)['start'])
        new_end = max(step_windows.tail(2)['end'])
        
        step_windows.drop(index=step_windows.index[-2:], inplace=True)
        
        delta_mean, p, r_con1, r_con2, r_between, con1_avg, con2_avg, avg_pval, con1_gini, con2_gini, gini_pval = window_calc(gene[gene['end'].between(new_start, new_end)])
        new_last_window = {'gene' : gene_name,
             'ID' : gene_annot.ID.unique()[0],
             'type' : gene_annot.type.unique()[0],
             'chr' : gene_annot.chr.unique()[0],
             'strand' : gene_annot.strand.unique()[0],          
             'start' : new_start,
             'end' : new_end,
             'x' : len(step_windows) + 1,
             'region': ','.join(gene[gene['end'].between(new_start, new_end)].region.unique()),
             'exon': ','.join(gene[gene['end'].between(new_start, new_end)].exon.astype(int).astype(str).unique()),
             'delta' : delta_mean, 'p' : p, 
             'r_con1' : r_con1, 'r_con2' : r_con2, 'r_between' : r_between,
             'con1_avg' : con1_avg, 'con2_avg' : con2_avg, 'avg_pval' : avg_pval, 
             'con1_gini' : con1_gini, 'con2_gini' : con2_gini, 'gini_pval' : gini_pval,
             'n': len(gene[gene['end'].between(new_start, new_end)])}
        
        new_last_window = pd.DataFrame(data=new_last_window, index = [0])
        step_windows = pd.concat([step_windows, new_last_window], ignore_index = True)
        
    #length calc
    pos_list = []
    for i in gene_annot.index:
        pos_list = pos_list + list(range(gene_annot.loc[i, 'start'] + 1, gene_annot.loc[i, 'end'] + 1))
    pos_list.sort()
    
    for i in step_windows.index:
        win_start = step_windows.loc[i, 'start']
        win_end = step_windows.loc[i, 'end']
        win_pos_list = [x for x in pos_list if x >= win_start and x<= win_end]
        step_windows.loc[i, 'len'] = len(win_pos_list)
    step_windows['len'] = step_windows.len.astype(int)

    step_output = pd.concat([step_output, step_windows], ignore_index = True)
    
    
    #SLIDE_WINDOWS: EACH WINDOW STARTS 10% OF WINDOW SIZE AFTER THE PREVIOUS WINDOW. OVERLAPPING WILL OCCUR, TAKEN CARE OF LATER.
    data = {}
    binSize = windowSize - 1 
    for i in range(0, len(gene), round(windowSize * 0.1) ):  
        delta_mean, p, r_con1, r_con2, r_between, con1_avg, con2_avg, avg_pval, con1_gini, con2_gini, gini_pval = window_calc(gene.loc[i:i+binSize, :])
        data[len(data)] = {'start' : min(gene.loc[i:i+binSize, :]['end']), #in gene df, 'end' coord is index1 inclusive position
             'end' : max(gene.loc[i:i+binSize, :]['end']),
             'region': ','.join(gene.loc[i:i+binSize, :].region.astype(str).unique()),
             'exon' : ','.join(gene.loc[i:i+binSize, :].exon.astype(int).astype(str).unique()),
             'delta' : delta_mean, 'p' : p, 
             'r_con1' : r_con1, 'r_con2' : r_con2, 'r_between' : r_between,
             'con1_avg' : con1_avg, 'con2_avg' : con2_avg, 'avg_pval' : avg_pval, 
             'con1_gini' : con1_gini, 'con2_gini' : con2_gini, 'gini_pval' : gini_pval,
             'n': len(gene.loc[i:i+binSize, :])}
              
    slide_windows = pd.DataFrame.from_dict(data=data, orient = 'index', columns=['gene', 'ID', 'type', 'chr', 'strand', 'region', 'exon', 'start', 'end',
                                                                           'x', 'delta', 'p', 'r_con1', 'r_con2', 'r_between', 
                                                                           'con1_avg', 'con2_avg', 'avg_pval', 
                                                                           'con1_gini', 'con2_gini', 'gini_pval', 'n', 'len'])
    slide_windows['gene'] = gene_name
    slide_windows['ID'] = gene_annot.ID.unique()[0]
    slide_windows['type'] = gene_annot.type.unique()[0]
    slide_windows['chr'] = gene_annot.chr.unique()[0] 
    slide_windows['strand'] = gene_annot.strand.unique()[0]
    slide_windows['err'] = ''

    #as opposed to the STEP calculation, where the last window was recombined with the prior window if it was too small, we will get rid of any very small widows at the end, as there will always be another window including those bases.
    slide_windows.drop(slide_windows[slide_windows['n'] < (0.5 * windowSize)].index, inplace=True)
              
 
    #length calc by creating a list of chromosomal coords in the canonical annotation (skipping introns) and then seeing how many of those positions are in each window. 
    pos_list = []
    for i in gene_annot.index:
        pos_list = pos_list + list(range(gene_annot.loc[i, 'start'] + 1, gene_annot.loc[i, 'end'] + 1))
    pos_list.sort()
    
    
    slide_windows_filtered = pd.DataFrame(columns=['gene', 'ID', 'type', 'chr', 'strand', 'region', 'exon', 'start', 'end',
                               'x', 'delta', 'p', 'r_con1', 'r_con2', 'r_between',
                               'con1_avg', 'con2_avg', 'avg_pval', 
                               'con1_gini', 'con2_gini', 'gini_pval', 'n', 'len'])
                               
    dedup_bool = [True for i in pos_list] #creates a bool list the length of the transcript
    for i, row in slide_windows.sort_values('delta', ascending = False).iterrows():
        if all(dedup_bool[pos_list.index(slide_windows.loc[i, 'start']): pos_list.index(slide_windows.loc[i, 'end'])]): #looks for overlapping windows in the deduplicated df. if all True, no positions overlap with an existing window
            slide_windows_filtered = pd.concat([slide_windows_filtered, row.to_frame().T], ignore_index = True, axis=0)
            for x in range(pos_list.index(slide_windows.loc[i, 'start']), pos_list.index(slide_windows.loc[i, 'end']) + 1):
                dedup_bool[x] = False #after appending the window, sets the corresponding values in the bool list to False to avoid overlaps with later windows
        
    for i in slide_windows_filtered.index:
        win_start = slide_windows_filtered.loc[i, 'start']
        win_end = slide_windows_filtered.loc[i, 'end']
        win_pos_list = [x for x in pos_list if x >= win_start and x<= win_end]
        slide_windows_filtered.loc[i, 'len'] = len(win_pos_list)
    slide_windows_filtered['len'] = slide_windows_filtered.len.astype(int)
    slide_windows_filtered['x'] = range(1, len(slide_windows_filtered) + 1)

    slide_output = pd.concat([slide_output, slide_windows_filtered], ignore_index = True)

os.chdir(outpath)
step_output.to_csv(str(windowSize) + 'nt_windows_step_' + chrom_in + '.txt', sep='\t', index=False, header = False)
step_output.to_csv(str(windowSize) + 'nt_windows_step_' + chrom_in + '.csv', index=False, header = False)
slide_output.to_csv(str(windowSize) + 'nt_windows_slide_' + chrom_in + '.txt', sep='\t', index=False, header = False)
slide_output.to_csv(str(windowSize) + 'nt_windows_slide_' + chrom_in + '.csv', index=False, header = False)
