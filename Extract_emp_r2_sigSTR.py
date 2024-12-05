#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 13:43:39 2024
Extract STR emperical R2 for genome-wide significant STR's
@author: pmerbaum
"""
import pandas as pd

sig = pd.read_csv('STR_meta/geno_sig_174k.txt', sep = '\t')
all_str = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/all_strata_all_chr_EHv500_ILL174K_b38_DF2_jan2024_meta_STR_sorted.txt.gz',sep = '\t') 

def transform_str(str_value, bp_value, t_value):
    return f'STR_ILL_174k_{bp_value}_T{t_value}'

def transform_C9_str(t_value):
    return f'STR_ILL_174k_27573528_T{t_value}'

sig_r2_combined = pd.DataFrame()
for chr in range(1, 23):
    filename = f'DF2_panel_EHv500_ILL174K_b38_DF2_jan2024_chr{chr}_empiricalR2_thresholded.txt'
    r2_info = pd.read_csv(filename, sep = '\t')

    bp_values = r2_info['STR'].str.extract(r'_([^_]+)_([^_]+)_([^_]+)_')
    r2_info['bp'] = bp_values[2]


# Apply the transformation to the 'STR' column of r2_chr{chr}
    r2_info['Threshold'] = r2_info['Threshold'].astype(int)
    special_str_mask = r2_info['STR'] == 'STR_ILL-174k_C9ORF72_GGCCCC'
    r2_info.loc[special_str_mask, 'STR_transformed'] = r2_info[special_str_mask].apply(lambda row: transform_C9_str(row['Threshold']), axis=1)
    r2_info.loc[~special_str_mask, 'STR_transformed'] = r2_info[~special_str_mask].apply(lambda row: transform_str(row['STR'], row['bp'], row['Threshold']), axis=1)

    merged_df = pd.merge(sig, r2_info, left_on='label', right_on='STR_transformed', how='inner')
    merged_df['chr'] = chr
    selected_columns = ['label','AF', 'bp_x','chr','r2', 'r2_spearman']
    sig_r2_combined = pd.concat([sig_r2_combined, merged_df[selected_columns]])
    
sig_r2_combined.reset_index(drop=True, inplace=True)
sig_r2_combined.to_csv('/hpc/hers_en/pmerbaum/str_imputation/post_meta_analysis/data/significance/emp_r2_from_NYGC_sig_STR.csv', sep='\t', index=False)

#weignted DR2 extraction

sig_r2_combined = pd.DataFrame()
def transform_str_no_t(str_value, bp_value):
    return f'STR_ILL_174k_{bp_value}'

def transform_C9_str_no_t():
    return 'STR_ILL_174k_27573528'

DR2_info = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/NYGC_ILL174K_weighted_DR2.txt', sep = '\t')
bp_values = DR2_info['chrom'].str.extract(r'_([^_]+)_([^_]+)_([^_]+)_')
DR2_info['bp']  = bp_values[2]
DR2_info['chr']  = bp_values[1]

special_str_mask_no_t = DR2_info['chrom'] == 'STR_ILL-174k_C9ORF72_GGCCCC'

DR2_info.loc[special_str_mask_no_t, 'STR_transformed'] = DR2_info[special_str_mask_no_t].apply(lambda row: transform_C9_str_no_t)
DR2_info.loc[~special_str_mask_no_t, 'STR_transformed'] = DR2_info[~special_str_mask_no_t].apply(lambda row: transform_str_no_t(row['chrom'], row['bp']), axis=1)


sig['bp'] = sig['label'].str.extract(r'_(\d+)_T')
sig['label_reformatted'] = 'STR_ILL_174k_' + sig['bp']
merged_df = pd.merge(sig, DR2_info, left_on='label_reformatted', right_on='STR_transformed', how='inner')
selected_columns = ['label','chr_x', 'Freq1', 'weighted_DR2', 'nonREF_AF']
sig_DR2 = pd.concat([sig_r2_combined, merged_df[selected_columns]])
sig_DR2.to_csv('sig_STR_NYGC_ILL174K_weighted_DR2.csv', index = False)

all_str['label_reformatted'] = 'STR_ILL_174k_' + all_str['POS'].astype(str)
merged_df = pd.merge(all_str, DR2_info, left_on='label_reformatted', right_on='STR_transformed', how='inner')
badly_imputed_thresholds = merged_df[merged_df['weighted_DR2'] <= 0.7]
badly_imputed = badly_imputed_thresholds.drop_duplicates(subset='label_reformatted')


columns_to_drop = ['Allele1', 'Allele2', 'chrom', 'start', 'end', 'RU', 'panel', 'mean', 'ID', 'bp','chr','STR_transformed', 'nalleles_common', 'hwep','var', 'MarkerName', 'Freq1', 'FreqSE', 'MinFreq', 'MaxFreq', 'Effect','StdErr', 'P-value', 'Direction', 'HetISq', 'HetChiSq', 'HetDf','HetPVal']
badly_imputed.drop(columns=columns_to_drop, inplace=True)
new_column_order = ['label_reformatted','CHR', 'POS','REF', 'ALT', 'weighted_DR2', 'r2', 'r2_spearman','nonREF_AF', 'nalleles']
badly_imputed = badly_imputed[new_column_order]
badly_imputed.rename(columns={'weighted_DR2': 'weighted_DR2_NYGC'}, inplace=True)
badly_imputed.rename(columns={'r2_spearman': 'r2_spearman_NYGC'}, inplace=True)
badly_imputed.rename(columns={'r2': 'r2_NYGC'}, inplace=True)
badly_imputed.to_csv('/Users/pmerbaum/Desktop/STR_weighted_DR2_below0.7_based_on_NYGC.txt',sep = '\t')

badly_imputedDF_thresholds = merged_df[merged_df['weighted_DR2'] <= 0.7]
badly_imputedDF = badly_imputedDF_thresholds.drop_duplicates(subset='label_reformatted')


columns_to_drop = ['Allele1', 'Allele2', 'chrom', 'start', 'end', 'RU', 'panel', 'mean', 'ID', 'bp','chr','STR_transformed', 'nalleles_common', 'hwep','var', 'MarkerName', 'Freq1', 'FreqSE', 'MinFreq', 'MaxFreq', 'Effect','StdErr', 'P-value', 'Direction', 'HetISq', 'HetChiSq', 'HetDf','HetPVal']
badly_imputedDF.drop(columns=columns_to_drop, inplace=True)
new_column_order = ['label_reformatted','CHR', 'POS','REF', 'ALT', 'weighted_DR2', 'r2', 'r2_spearman','nonREF_AF', 'nalleles']
badly_imputedDF = badly_imputedDF[new_column_order]
badly_imputedDF.rename(columns={'weighted_DR2': 'weighted_DR2_DF2'}, inplace=True)
badly_imputedDF.rename(columns={'r2_spearman': 'r2_spearman_DF2'}, inplace=True)
badly_imputedDF.rename(columns={'r2': 'r2_DF2'}, inplace=True)

weight = pd.merge(badly_imputed, badly_imputedDF, on = ['label_reformatted','CHR', 'POS','REF', 'ALT', 'nonREF_AF', 'nalleles'])
#STR's that are in my meta analysis

#all badly imputed STRs
NYGC_info = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/NYGC_ILL174K_weighted_DR2.txt', sep = '\t')
DF2_info = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/DF2_ILL174K_weighted_DR2.txt', sep = '\t')
merged = pd.merge(DF2_info, NYGC_info, on = ['chrom','start', 'end', 'RU', 'ID', 'panel', 'REF', 'ALT','nonREF_AF', 'nalleles', 'nalleles_common', 'hwep', 'mean', 'var'])
merged.drop(columns=['hwep', 'mean', 'var', 'nalleles_common', 'start', 'end', 'RU', 'ID'], inplace=True)
merged.rename(columns={'weighted_DR2_x': 'weighted_DR2_DF2'}, inplace=True)
merged.rename(columns={'weighted_DR2_y': 'weighted_DR2_NYGC'}, inplace=True)
merged.rename(columns={'r2_x': 'r2_DF2'}, inplace=True)
merged.rename(columns={'r2_y': 'r2_NYGC'}, inplace=True)
merged.rename(columns={'r2_spearman_x': 'r2_spearman_DF2'}, inplace=True)
merged.rename(columns={'r2_spearman_y': 'r2_spearman_NYGC'}, inplace=True)
merged.rename(columns={'chrom': 'STR_ID'}, inplace=True)
badly_imputed_merged = merged[merged['weighted_DR2_DF2'] <= 0.7]
badly_imputed_merged.to_csv('/Users/pmerbaum/Desktop/STR_weighted_DR2_below0.7.txt',sep = '\t')

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

emp_df2 = pd.read_csv('emp_r2_from_NYGC_sig_STR.csv', sep ='\t')
# Assuming 'emp_df2' is your DataFrame
plt.scatter(emp_df2['chr'], emp_df2['r2_spearman'])
plt.xlabel('STR')
plt.ylabel('empirical R2')
plt.title('STR empirical R2')
plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))

# Add horizontal line at y = 0.7
plt.axhline(y=0.7, color='r', linestyle='--', label='Cut-off')
plt.legend()  # Show legend with the cut-off line

plt.xticks(rotation=90)  # Rotate x-axis labels for better readability
plt.tight_layout()  # Adjust layout to prevent clipping of labels
plt.show()

