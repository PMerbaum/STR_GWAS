#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 12:59:22 2024
Maiami plot 

@author: pmerbaum
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

    

base5 = pd.read_csv('/Users/pmerbaum/Desktop/Phd_docs/Base5/s.meta.hg19.txt', sep = '\t')
gwas2018_euro_only = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/ALS_sumstats_EUR_only.txt.gz', sep = '\t')

#gwas_subset = gwas2018_euro_only[(gwas2018_euro_only['chr'] == 6) & 
#                                   (gwas2018_euro_only['bp'] >= 28477063) & 
#                                  (gwas2018_euro_only['bp'] <= 33448936)] 
gwas_subset = gwas2018_euro_only[(gwas2018_euro_only['chr'] == 6) & 
                                   (gwas2018_euro_only['bp'] >= 32000000) & 
                                   (gwas2018_euro_only['bp'] <= 33000000)] #to correspond with base5 data

base5_subset = base5[(base5['POS_hg19'] >= 32000000) & 
                                   (base5['POS_hg19'] <= 33000000)]

base5_subset['-log10P'] = -np.log10(base5['P.value'])
gwas_subset['-log10P'] = -np.log10(gwas_subset['P-value'])

threshold = 5e-08
threshold_logp = -np.log10(threshold)
base5_subset['Color'] = np.where(base5_subset['P.value'] < threshold, 'coral', 'lightblue')
gwas_subset['Color'] = np.where(gwas_subset['P-value'] < threshold, 'coral', 'lightblue')


fig, ax = plt.subplots(figsize=(10, 8))

# Plot meta-analysis
ax.scatter(base5_subset['POS_hg19'], base5_subset['-log10P'], c=base5_subset['Color'], cmap='tab20', alpha=0.6, edgecolors="w", s=50)
ax.set_ylabel('-log10(P) (BASE5 6 strata meta-analysis)', fontsize=12)
ax.axhline(y=-threshold_logp, color='gray', linestyle='--', alpha=0.7, label=f'Significance threshold: P < {threshold}')
#ax.set_title('Miami Plot', fontsize=16)


# Plot GWAS results
ax.scatter(gwas_subset['bp'], -gwas_subset['-log10P'], c=gwas_subset['Color'], cmap='tab20', alpha=0.6, edgecolors="w", s=50)
ax.set_ylabel('-log10(P) (ALS GWAS 2021)', fontsize=12)
ax.axhline(y=threshold_logp, color='gray', linestyle='--', alpha=0.7, label=f'Significance threshold: P < {threshold}')
ax.set_xlabel('Base Pair Position', fontsize=12)
ax.grid(False)

# Adjust layout and show plot
plt.tight_layout()
plt.show()

##############################################################################
#p-p plot between base5 and gwas results
gwas_subset['bp'] = gwas_subset['bp'].astype(int)
only_eur = pd.merge(base5, gwas_subset, left_on='POS_hg19', right_on='bp', how='left')

plt.figure(figsize=(8, 6))
sns.scatterplot(x='-log10P_x', y='-log10P_y', data=only_eur)

# Add a trendline
max_value = max(only_eur['-log10P_x'].max(), only_eur['-log10P_y'].max())
plt.plot([0, max_value], [0, max_value], color='gray', linestyle='--')

# Mark significant p-values only in one dataset with red dots
significant_only_in_x = (only_eur['-log10P_x'] < -np.log10(5e-08)) & ~(only_eur['-log10P_y'] < -np.log10(5e-08))
significant_only_in_y = (only_eur['-log10P_y'] < -np.log10(5e-08)) & ~(only_eur['-log10P_x'] < -np.log10(5e-08))

# Plot significant points
plt.scatter(only_eur.loc[significant_only_in_x, '-log10P_x'], 
            only_eur.loc[significant_only_in_x, '-log10P_y'], 
            color='red', label='Significant only in BASE5 meta')

plt.scatter(only_eur.loc[significant_only_in_y, '-log10P_x'], 
            only_eur.loc[significant_only_in_y, '-log10P_y'], 
            color='pink', label='Significant only in GWAS 2018')

# Add labels and title
plt.xlabel('-logP BASE5 meta')
plt.ylabel('-logP GWAS 2018')
plt.title('P-p plot for BASE5 meta and GWAS 2018 EUR only')

significant_x_marker_ids = only_eur.loc[significant_only_in_x, 'MarkerName_x']
significant_y_marker_ids = only_eur.loc[significant_only_in_y, 'MarkerName_y']