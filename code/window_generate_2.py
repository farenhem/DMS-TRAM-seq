#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
import math
import datetime
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--size", dest="SIZE", help="Integer-value size of windows to be calculated.")
parser.add_option("--path", dest="PATH", help="Full path to working directory, which contains the bedGraph and set_windows folders.")
parser.add_option("--chr", dest="CHR", help="Chromosome name. Ex: chr1")
parser.add_option("--comparison", dest="COMP", help="The prefixes of both sample conditions, connected by an underscore. Ex: condition1_condition2.")

    
(options, args) = parser.parse_args()

wkdir = str(options.PATH)
windowSize = int(options.SIZE)
chrom = str(options.CHR)
comp = str(options.COMP)

print('Starting annotation for ' + chrom + '.')

os.chdir(wkdir + '/set_windows/chrom')

STEP_IDs_pos = pd.read_csv(comp + '_' + chrom + '_' + str(windowSize) + 'bp_deltas_STEP_annot_pos.txt', sep ='\t', header=None, names = ['chr', 'start', 'end', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval', 'ID', 'exon', 'HepG2', 'K562'])
STEP_IDs_neg = pd.read_csv(comp + '_' + chrom + '_' + str(windowSize) + 'bp_deltas_STEP_annot_neg.txt', sep ='\t', header=None, names = ['chr', 'start', 'end', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval', 'ID', 'exon', 'HepG2', 'K562'])
STEP_IDs_pos['strand'] = '+'
STEP_IDs_neg['strand'] = '-'
STEP = (pd.concat([STEP_IDs_pos, STEP_IDs_neg])).reset_index(drop=True)

SLIDE_IDs_pos = pd.read_csv(comp + '_' + chrom + '_' +str(windowSize) + 'bp_deltas_SLIDE_annot_pos.txt', sep ='\t', header=None, names = ['chr', 'start', 'end', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval', 'ID', 'exon', 'HepG2', 'K562'])
SLIDE_IDs_neg = pd.read_csv(comp + '_' + chrom + '_' +str(windowSize) + 'bp_deltas_SLIDE_annot_neg.txt', sep ='\t', header=None, names = ['chr', 'start', 'end', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval', 'ID', 'exon', 'HepG2', 'K562'])
SLIDE_IDs_pos['strand'] = '+'
SLIDE_IDs_neg['strand'] = '-'
SLIDE = (pd.concat([SLIDE_IDs_pos, SLIDE_IDs_neg])).reset_index(drop=True)


os.chdir(wkdir+'/reference_annot')
genedf = (pd.read_csv('ensemblID_genename.txt', sep ='\t', names=['ID', 'gene'])).fillna('')
ID_type = (pd.read_csv('ensemblID_transcripttype.txt', sep='\t', names=['ID', 'type'])).fillna('')

STEP['gene']  = ''
STEP['type'] = ''
for i in range(len(STEP)):
    name_list = []
    type_list = []
    ID_list = str(STEP['ID'][i]).split(';')
    if len(ID_list) > 2:
        continue
    if ID_list != ['']:
        for x in ID_list:
            if x in list(genedf['ID']):
                name = genedf.loc[genedf['ID'] == x, 'gene'].item()
                if name != '':
                        name_list.append(str(name))
                ttype = ID_type.loc[ID_type['ID'] == x, 'type'].item()
                if ttype != '':
                        type_list.append(str(ttype))
        name_string = ';'.join(name_list)
        type_string = ';'.join(type_list)
        STEP.loc[i, 'gene'] = name_string
        STEP.loc[i, 'type'] = type_string
    
SLIDE['gene']  = ''
SLIDE['type'] = ''
for i in range(len(SLIDE)):
    name_list = []
    type_list = []
    ID_list = str(SLIDE['ID'][i]).split(';')
    if len(ID_list) > 2:
        continue
    if ID_list != ['']:
        for x in ID_list:
            if x in list(genedf['ID']): #needed in case annotations don't exactly match (rare but will throw an error with .item function below)
                name = genedf.loc[genedf['ID'] == x, 'gene'].item()
                if name != '':
                        name_list.append(str(name))
                ttype = ID_type.loc[ID_type['ID'] == x, 'type'].item()
                if ttype != '':
                        type_list.append(str(ttype))
        name_string = ';'.join(name_list)
        type_string = ';'.join(type_list)
        SLIDE.loc[i, 'gene'] = name_string
        SLIDE.loc[i, 'type'] = type_string


os.chdir(wkdir+'/reference_annot')
SG_data = pd.read_csv('SG_data_simple.csv', names = ['ID', 'SG_FC', 'SG_localization'])                                                                                                         
SLIDE = SLIDE.merge(SG_data, on = 'ID', how ='left')
STEP = STEP.merge(SG_data, on = 'ID', how ='left')        

os.chdir(wkdir+'/set_windows/chrom')
STEP.to_csv(comp + '_' + chrom + '_' + str(windowSize) + 'bp_STEP_deltas.csv', index = False, header = False, columns = ['chr', 'start', 'end', 'ID', 'gene', 'type', 'exon', 'strand', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval', 'HepG2', 'K562', 'SG_FC', 'SG_localization'])
SLIDE.to_csv(comp + '_' + chrom + '_' + str(windowSize) + 'bp_SLIDE_deltas.csv', index = False, header = False, columns = ['chr', 'start', 'end', 'ID', 'gene', 'type', 'exon', 'strand', 'delta', 'p', 'con1_avg', 'con2_avg', 'avg_pval', 'con1_gini', 'con2_gini', 'gini_pval', 'HepG2', 'K562', 'SG_FC', 'SG_localization'])

print('Finished annotation for ' + chrom + '.')