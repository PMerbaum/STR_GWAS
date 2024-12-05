#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 11:10:00 2024

@author: pmerbaum
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

data = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/geno_sig_174k.txt', sep = '\t')

info = [(row['chr'], row['bp']) for _, row in data.iterrows()] #tuple containing chr and bp 
info_set = set(info) #unique values 

meta = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/all_strata_all_chr_EHv500_ILL174K_b38_DF2_meta_exclusion.txt', sep = '\t')

chr_values = [tup[0] for tup in info_set]
pos_values = [tup[1] for tup in info_set]

# Create a boolean mask to filter rows where CHR and POS match
mask = (meta['CHR'].isin(chr_values)) & (meta['POS'].isin(pos_values))
filtered_marker_names = meta.loc[mask, 'MarkerName'] #all STR thresholds for genome-wide significant STR's
mask2 = meta['MarkerName'].isin(filtered_marker_names)
all_th_sigSTR = meta[mask2] #df with all info
all_th_sigSTR['STR'].unique() #these are my STR's of interest

split = all_th_sigSTR['MarkerName'].str.rsplit('_', n=1, expand=True)
all_th_sigSTR['STR'] = split[0]
all_th_sigSTR['Threshold'] = split[1]
all_th_sigSTR['Threshold'] = all_th_sigSTR['Threshold'].str.extract(r'T(\d+)').astype(int)
#all_th_sigSTR['lower_ci'] = all_th_sigSTR['Effect'] - 1.96 * all_th_sigSTR['StdErr']
#all_th_sigSTR['upper_ci'] = all_th_sigSTR['Effect'] + 1.96 * all_th_sigSTR['StdErr']


#plot each STR
STR = all_th_sigSTR[all_th_sigSTR['STR'] == 'STR_ILL_174k_49903902']
colors = np.where(STR['P-value'] < 5e-08, 'orange', 'blue')

plt.errorbar(STR['Threshold'], np.exp(STR['Effect']), yerr=np.exp(STR['Effect']) * (np.exp(STR['StdErr']) - 1), fmt='o', label='Effect', zorder=1)
plt.scatter(STR['Threshold'], np.exp(STR['Effect']), label='Effect', color=colors, zorder=2)

legend_handles = [mpatches.Patch(color='blue', label='Non-genome-wide significant'), mpatches.Patch(color='orange', label='Genome-wide significant')]
# Set plot labels and title
plt.xlabel('Threshold')
plt.ylabel('log Odds Ratio ')
plt.title('STR_ILL_174k_49903902 log(OR) per threshold')

plt.legend(handles=legend_handles)
# Show plot
plt.show()
