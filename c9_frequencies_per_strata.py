#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 11:27:18 2024

@author: pmerbaum
"""
import pandas as pd
import numpy as np
#s = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/s1.no_overlap.hm3_panel_EHv500_ILL174K_b38_DF2_jan2024_chr9_saige_step2_assoc_recoded.txt' , sep = '\s+')

s_cohorts = [1,2,3,4,6,7]

c9_lengths = pd.DataFrame()
for s in s_cohorts:
    filename = f's{s}.no_overlap.hm3_panel_EHv500_ILL174K_b38_DF2_jan2024_chr9_saige_step2_assoc_recoded.txt'
    s_info = pd.read_csv(filename, sep = '\s+')
    STR_info = s_info['SNPID'].str.extract(r'_([^_]+)_([^_]+)_([^_]+)_([^_]+)')
    
    df = pd.DataFrame({'threshold': STR_info.iloc[:, 3]})
    df['threshold'] = df['threshold'].str.extract(r'T(\d+)').astype(int)
    df['AF'] = s_info['AF_Allele2']
    df['STR'] = s_info['SNPID']
    df['strata'] = s
    
    c9_lengths = pd.concat([c9_lengths, df])
    

c9_only = c9_lengths[c9_lengths['STR'].str.match(r'STR_ILL_174k_27573528_T\d+')]

import matplotlib.pyplot as plt

y_ticks = np.arange(0, 1.1, 0.1)

plt.figure(figsize=(10, 6))
for strata, color in zip(c9_only['strata'].unique(), ['blue', 'green', 'red', 'orange', 'purple', 'brown']):
    plt.scatter(c9_only[c9_only['strata'] == strata]['threshold'],
                c9_only[c9_only['strata'] == strata]['AF'],
                color=color,
                label=f'Strata {strata}')

plt.ylim(0, 1)
plt.yticks(y_ticks)
# Add labels and title
plt.xlabel('Threshold')
plt.ylabel('Allele frequency')
plt.title('C9 thresholds and allele frequencies per strata')
plt.legend()

# Show plot
plt.grid(True)
plt.tight_layout()
plt.show()