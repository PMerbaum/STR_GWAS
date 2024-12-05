#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 11:43:15 2024

@author: pmerbaum
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sig = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/geno_sig_174k.txt') #filtered_data
maf_data = pd.DataFrame(columns=['label'])

s_number = [1, 2, 3, 4, 6, 7]
for s in s_number:
    strata = pd.read_csv(f's{s}_all_chr_EHv500_ILL174K_b38_DF2_format.txt', sep = ',')
    strata = strata.drop_duplicates('SNPID')
    merged_data = pd.merge(sig, strata[['SNPID', 'AF_Allele2']], how='inner', left_on='label', right_on='SNPID')
    merged_data.rename(columns={'AF_Allele2': f'{s}_MAF'}, inplace=True)
    maf_data = pd.merge(maf_data, merged_data[['label', f'{s}_MAF']], on='label', how='outer')
    
    

maf_data.to_csv('maf_data_per_strata.txt', sep='\t', index=False)