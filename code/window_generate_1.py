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
parser.add_option("--size", dest="SIZE", help="Integer-value size of windows to be calculated.")
parser.add_option("--path", dest="PATH", help="Full path to working directory, which contains the bedGraph and set_windows folders.")
parser.add_option("--condition1", dest="CON1", help="Prefix to the sample names of the first condition.")
parser.add_option("--condition2", dest="CON2", help="Prefix to the sample names of the second condition.")
    
(options, args) = parser.parse_args()

wkdir = str(options.PATH)
windowSize = int(options.SIZE)
binSize = windowSize - 1
chrom = str(options.CHR)
con1 = str(options.CON1)
con2 = str(options.CON2)

print('Starting window calculation for ' + chrom + '.')

sample1 = con1 + "1"
sample2 = con1 + "2"
sample3 = con1 + "3"
sample4 = con2 + "1"
sample5 = con2 + "2"
sample6 = con2 + "3"


os.chdir(wkdir+'/bedGraph/chrom')

sample1_neg = pd.read_csv(sample1 + '_' + chrom + '_mmRate_neg.bedGraph', sep='\t', header = None, names = ['chr', 'start', 'pos', 'mm'])
sample2_neg = pd.read_csv(sample2 + '_' + chrom + '_mmRate_neg.bedGraph', sep='\t', header = None, names = ['chr', 'start', 'pos', 'mm'])
sample3_neg = pd.read_csv(sample3 + '_' + chrom + '_mmRate_neg.bedGraph', sep='\t', header = None, names = ['chr', 'start', 'pos', 'mm'])
sample4_neg = pd.read_csv(sample4 + '_' + chrom + '_mmRate_neg.bedGraph', sep='\t', header = None, names = ['chr', 'start', 'pos', 'mm'])
sample5_neg = pd.read_csv(sample5 + '_' + chrom + '_mmRate_neg.bedGraph', sep='\t', header = None, names = ['chr', 'start', 'pos', 'mm'])
sample6_neg = pd.read_csv(sample6 + '_' + chrom + '_mmRate_neg.bedGraph', sep='\t', header = None, names = ['chr', 'start', 'pos', 'mm'])


sample1_pos = pd.read_csv(sample1 + '_' + chrom + '_mmRate_pos.bedGraph', sep='\t', header = None, names = ['chr', 'start', 'pos', 'mm'])
sample2_pos = pd.read_csv(sample2 + '_' + chrom + '_mmRate_pos.bedGraph', sep='\t', header = None, names = ['chr', 'start', 'pos', 'mm'])
sample3_pos = pd.read_csv(sample3 + '_' + chrom + '_mmRate_pos.bedGraph', sep='\t', header = None, names = ['chr', 'start', 'pos', 'mm'])
sample4_pos = pd.read_csv(sample4 + '_' + chrom + '_mmRate_pos.bedGraph', sep='\t', header = None, names = ['chr', 'start', 'pos', 'mm'])
sample5_pos = pd.read_csv(sample5 + '_' + chrom + '_mmRate_pos.bedGraph', sep='\t', header = None, names = ['chr', 'start', 'pos', 'mm'])
sample6_pos = pd.read_csv(sample6 + '_' + chrom + '_mmRate_pos.bedGraph', sep='\t', header = None, names = ['chr', 'start', 'pos', 'mm'])


pos = pd.DataFrame()
pos['chr'] = sample4_pos['chr']
pos['pos'] = sample4_pos['pos']

pos['sample4'] = sample4_pos['mm']
pos['sample5'] = sample5_pos['mm']
pos['sample6'] = sample6_pos['mm']
pos['sample1'] = sample1_pos['mm']
pos['sample2'] = sample2_pos['mm']
pos['sample3'] = sample3_pos['mm']


neg = pd.DataFrame()
neg['chr'] = sample4_neg['chr']
neg['pos'] = sample4_neg['pos']

neg['sample4'] = sample4_neg['mm']
neg['sample5'] = sample5_neg['mm']
neg['sample6'] = sample6_neg['mm']
neg['sample1'] = sample1_neg['mm']
neg['sample2'] = sample2_neg['mm']
neg['sample3'] = sample3_neg['mm']


pos.fillna('', inplace = True)
neg.fillna('', inplace = True)

def gini(x):
    #calculates gini coeff for an array/list/column
    mg = np.mean(x)
    if (mg == 0):
      g = 0
    else:
      mad = np.abs(np.subtract.outer(x, x)).mean()
      g = 0.5 * mad/np.mean(x)
    return g
    
def BH_pval_correction(pvals_in):
    pvals = np.array(pvals_in)
    ranks = rankdata(pvals, method='ordinal')
    padj = np.empty(len(pvals))
    p_m_over_i = pvals * len(pvals) / ranks
    for i in range(0, len(ranks)):
        current_rank = ranks[i]
        padj[i] = min(1, min(p_m_over_i[ranks >= current_rank])) #ensures no p-val of a higher rank is lower; max of 1.
    return padj 

if windowSize >= 40: #define upper winsorization limit based on window size to ensure removal of at least 1 potential outlier
    upper = 0.05 #5% upper wins limit is the typical value
else:
    upper = 0.10 #10% for small windows; less than 20 bp may still be problematic!

def get_deltas(df): #input df will already be sliced to only the window of interest
    #find delta values between the two conditions, paired by replicate
    if len(df) < windowSize: #skips small windows (occurs at end of chrom), prevents error in pearsonr function for those with len < 2. these rows will be deleted later via dropna. 
        delta_mean = np.nan
        p = np.nan
        con1_avg = np.nan
        con2_avg = np.nan
        avg_pval = np.nan
        con1_gini = np.nan
        con2_gini = np.nan
        gini_pval = np.nan
    elif ((df.sample1.values >= 0.01).sum() / windowSize) < 0.3: #at least 30% of datapoints must have mmRate > 0.01; this will skip incorrectly-stranded regions
        delta_mean = np.nan
        p = np.nan
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
        
        rbetween.append(pearsonr(df.sample1.values, df.sample4.values).statistic)
        rbetween.append(pearsonr(df.sample1.values, df.sample5.values).statistic)
        rbetween.append(pearsonr(df.sample1.values, df.sample6.values).statistic)
        rbetween.append(pearsonr(df.sample2.values, df.sample4.values).statistic)
        rbetween.append(pearsonr(df.sample2.values, df.sample5.values).statistic)
        rbetween.append(pearsonr(df.sample2.values, df.sample6.values).statistic)
        rbetween.append(pearsonr(df.sample3.values, df.sample4.values).statistic)
        rbetween.append(pearsonr(df.sample3.values, df.sample5.values).statistic)
        rbetween.append(pearsonr(df.sample3.values, df.sample6.values).statistic)

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
        
        sample1_gini = gini(stats.mstats.winsorize(df.sample1.values, limits=[0, upper]))
        sample2_gini = gini(stats.mstats.winsorize(df.sample2.values, limits=[0, upper]))
        sample3_gini = gini(stats.mstats.winsorize(df.sample3.values, limits=[0, upper]))
        sample4_gini = gini(stats.mstats.winsorize(df.sample4.values, limits=[0, upper]))
        sample5_gini = gini(stats.mstats.winsorize(df.sample5.values, limits=[0, upper]))
        sample6_gini = gini(stats.mstats.winsorize(df.sample6.values, limits=[0, upper]))
        
        con1_gini = (sample1_gini + sample2_gini + sample3_gini) / 3
        con2_gini = (sample4_gini + sample5_gini + sample6_gini) / 3
        gini_pval = ttest_ind([sample1_gini, sample2_gini, sample3_gini], [sample4_gini, sample5_gini, sample6_gini]).pvalue
        
    return delta_mean, p, con1_avg, con2_avg, avg_pval, con1_gini, con2_gini, gini_pval



neg_data = {}
for i in range(0, len(neg), windowSize): #next start coord is from the row after the prior window ends
    delta, p, con1_avg, con2_avg, avg_pval, con1_gini, con2_gini, gini_pval = get_deltas(neg.loc[i:i+binSize, :])
    neg_data[len(neg_data)] = {'chr' : neg.at[i, 'chr'],
         'start_coord' : neg.at[i, 'pos'],
         'end_coord' : neg.at[min(i+binSize, len(neg)-1), 'pos'],
         'delta' : delta,
         'p' : p, 
         'con1_avg' : con1_avg,
         'con2_avg'  : con2_avg,
         'avg_pval' : avg_pval,
         'con1_gini' : con1_gini,
         'con2_gini' : con2_gini,
         'gini_pval' : gini_pval}
            
neg_windows = pd.DataFrame.from_dict(data=neg_data, orient='index')


pos_data = {}
for i in range(0, len(pos), windowSize): 
    delta, p, con1_avg, con2_avg, avg_pval, con1_gini, con2_gini, gini_pval  = get_deltas(pos.loc[i:i+binSize, :])
    pos_data[len(pos_data)] = {'chr' : pos.at[i, 'chr'],
         'start_coord' : pos.at[i, 'pos'],
         'end_coord' : pos.at[min(i+binSize, len(pos)-1), 'pos'],
         'delta' : delta,
         'p' : p, 
         'con1_avg' : con1_avg,
         'con2_avg'  : con2_avg,
         'avg_pval' : avg_pval,
         'con1_gini' : con1_gini,
         'con2_gini' : con2_gini,
         'gini_pval' : gini_pval}
            
pos_windows = pd.DataFrame.from_dict(data=pos_data, orient='index')




#remove windows without deltas (small windows at end of chr) 
neg_windows.dropna(axis=0, subset = ['delta'], inplace=True)
neg_windows.reset_index(drop=True, inplace=True)
pos_windows.dropna(axis=0, subset = ['delta'], inplace=True)
pos_windows.reset_index(drop=True, inplace=True)


os.chdir(wkdir+'/set_windows/chrom')

neg_windows.to_csv(con1 + '_' + con2 + '_' + chrom + '_' + str(windowSize) + 'bp_deltas_STEP_neg.txt', sep='\t', index = False, header = False, columns = ['chr', 'start_coord', 'end_coord', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval'])
pos_windows.to_csv(con1 + '_' + con2 + '_' + chrom + '_' + str(windowSize) + 'bp_deltas_STEP_pos.txt', sep='\t', index = False, header = False, columns = ['chr', 'start_coord', 'end_coord', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval'])



neg_data = {}
for i in range(0, len(neg)): #SLIDING/ROLLING- START COORD ONLY DIFFERS BY 1 BP
    delta, p, con1_avg, con2_avg, avg_pval, con1_gini, con2_gini, gini_pval  = get_deltas(neg.loc[i:i+binSize, :])
    neg_data[len(neg_data)] = {'chr' : neg.at[i, 'chr'],
         'start_coord' : neg.at[i, 'pos'],
         'end_coord' : neg.at[min(i+binSize, len(neg)-1), 'pos'],
         'delta' : delta,
         'p': p, 
         'con1_avg' : con1_avg,
         'con2_avg'  : con2_avg,
         'avg_pval' : avg_pval,
         'con1_gini' : con1_gini,
         'con2_gini' : con2_gini,
         'gini_pval' : gini_pval}
            
neg_windows_sliding = pd.DataFrame.from_dict(data=neg_data, orient='index')
neg_windows_sliding.dropna(axis=0, subset = ['delta'], inplace=True)
neg_windows_sliding.reset_index(drop=True, inplace=True)#remove small windows



pos_data = {}
for i in range(0, len(pos)): 
    delta, p, con1_avg, con2_avg, avg_pval, con1_gini, con2_gini, gini_pval  = get_deltas(pos.loc[i:i+binSize, :])
    pos_data[len(pos_data)] = {'chr' : pos.at[i, 'chr'],
         'start_coord' : pos.at[i, 'pos'],
         'end_coord' : pos.at[min(i+binSize, len(pos)-1), 'pos'],
         'delta' : delta,
         'p' : p, 
         'con1_avg' : con1_avg,
         'con2_avg'  : con2_avg,
         'avg_pval' : avg_pval,
         'con1_gini' : con1_gini,
         'con2_gini' : con2_gini,
         'gini_pval' : gini_pval}
            
pos_windows_sliding = pd.DataFrame.from_dict(data=pos_data, orient='index')
pos_windows_sliding.dropna(axis=0, inplace=True, subset=['delta'])
pos_windows_sliding.reset_index(drop=True, inplace=True) #remove small windows


neg_slide_dedup = pd.DataFrame(columns = ['chr', 'start_coord', 'end_coord', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval'])
pos_slide_dedup = pd.DataFrame(columns = ['chr', 'start_coord', 'end_coord', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval'])


dedup_bool = [True for i in range(max(neg_windows_sliding['end_coord']))] #creates a bool list the length of the chromosome, all values start as True but will become False when a window covering that region is included in the final dataframe
for i, row in neg_windows_sliding.sort_values('delta', ascending = False).iterrows():
    if all(dedup_bool[neg_windows_sliding.loc[i, 'start_coord'] - 1 : neg_windows_sliding.loc[i, 'end_coord']]): #looks for overlapping windows in the deduplicated df. if True, no positions overlap with an existing window
        neg_slide_dedup = pd.concat([neg_slide_dedup, row.to_frame().T], ignore_index = True, axis=0)
        for x in range(neg_windows_sliding.loc[i, 'start_coord'] - 1 , neg_windows_sliding.loc[i, 'end_coord']):
            dedup_bool[x] = False #after appending the window, sets the corresponding values in the bool list to False to avoid overlaps with later windows


dedup_bool = [True for i in range(max(pos_windows_sliding['end_coord']))] #creates a bool list the length of the chromosome, all values start as True but will become False when a window covering that region is included in the final dataframe
for i, row in pos_windows_sliding.sort_values('delta', ascending = False).iterrows():
    if all(dedup_bool[pos_windows_sliding.loc[i, 'start_coord'] - 1 : pos_windows_sliding.loc[i, 'end_coord']]): #looks for overlapping windows in the deduplicated df. if True, no positions overlap with an existing window
        pos_slide_dedup = pd.concat([pos_slide_dedup, row.to_frame().T], ignore_index = True, axis=0)
        for x in range(pos_windows_sliding.loc[i, 'start_coord'] - 1 , pos_windows_sliding.loc[i, 'end_coord']):
            dedup_bool[x] = False #after appending the window, sets the corresponding values in the bool list to False to avoid overlaps with later windows       

neg_slide_dedup.to_csv(con1 + '_' + con2 + '_' + chrom + '_' + str(windowSize) + 'bp_deltas_SLIDE_neg.txt', sep='\t', index = False, header = False, columns = ['chr', 'start_coord', 'end_coord', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval'])
pos_slide_dedup.to_csv(con1 + '_' + con2 + '_' + chrom + '_' + str(windowSize) + 'bp_deltas_SLIDE_pos.txt', sep='\t', index = False, header = False, columns = ['chr', 'start_coord', 'end_coord', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval'])

print('Finished window calculation for ' + chrom + '.')