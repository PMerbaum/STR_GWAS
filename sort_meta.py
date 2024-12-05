#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:07:24 2024

@author: pmerbaum
sort files for Locus Zoom

"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def transform_marker_id(marker_id):
    parts = marker_id.split(':')
    chromosome = 'chr' + parts[0]
    return chromosome + ':' + parts[1]

bad_snps = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/SNP_genotype_DR2_filtered_0.7_MAF.txt.gz', sep = '\t')
bad_snps['SNP_ID'] = bad_snps['ID'].str.split(':').str[:2].str.join(':')


unsorted_SNP = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/all_strata_all_chr_genotype_EHv500_ILL174K_b38_DF2_jan2024_meta_SNPs.txt.gz', sep = '\t')
unsorted_STR = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/all_strata_all_chr_genotype_EHv500_ILL174K_b38_DF2_jan2024_meta_STRs.txt.gz', sep = '\t')
#sig = unsorted_STR[unsorted_STR['P-value'] < 5e-08]
#unique = sig[['POS', 'CHR']].drop_duplicates()
#sig.to_csv('/Users/pmerbaum/Desktop/STR_meta/genome_wide_sig_STR_genotype.txt.gz', sep='\t', compression='gzip', index=False)


#change to formats chr:pos or chr:pos_ref/alt
unsorted_SNP['CHR'] = unsorted_SNP['MarkerName'].str.extract(r'chr(\d+)')
unsorted_SNP['POS'] = unsorted_SNP['MarkerName'].str.extract(r'chr\d+_(\d+)')
unsorted_SNP['POS'] = unsorted_SNP['POS'].astype(int)
unsorted_SNP['Marker_ID'] = unsorted_SNP['CHR'].astype(str) + ':' + unsorted_SNP['POS'].astype(str)
unsorted_SNP['SNP_ID'] = unsorted_SNP['Marker_ID'].apply(transform_marker_id)
unsorted_SNP= unsorted_SNP[~unsorted_SNP['SNP_ID'].isin(bad_snps['SNP_ID'])]


unsorted_STR['threshold'] = unsorted_STR['MarkerName'].str.extract(r'T(\d+)')
unsorted_STR['Marker_ID'] = unsorted_STR['CHR'].astype(str) + ':' + unsorted_STR['POS'].astype(str) + '_s/T' + unsorted_STR['threshold'].astype(str)

#exclude STRs that failed QC
exclude = [27302982, 27999436, 52861004, 52886977,3850273]
unsorted_STR_filtered = unsorted_STR[~unsorted_STR['POS'].isin(exclude)]

joined_STR_SNP_unsorted = pd.concat([unsorted_STR_filtered, unsorted_SNP],ignore_index=True)

joined_STR_SNP_sorted = joined_STR_SNP_unsorted.sort_values(by=['CHR', 'POS'])
joined_STR_SNP_sorted.to_csv('/Users/pmerbaum/Desktop/STR_meta/all_strata_all_chr_genotype_EHv500_ILL174K_b38_DF2_jan2024_meta_STR_SNP_nobad_SNP_nobadSTR_sorted.txt.gz', sep='\t', compression='gzip', index=False)

#############################################################################
rsID_hg38 = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/STR_GWAS_SNP_RSids.txt.gz')
rsID_hg38['RS'] = 'rs' + rsID_hg38['RS'].astype(str)
sorted_SNP_rs = pd.merge(unsorted_SNP, rsID_hg38, left_on='Marker_ID', right_on='Chrom_Pos', how='left') #contains duplicates from different versions of dbSNP.hg38. Keep the older snp by executing next line 
sorted_SNP_rs_unique = sorted_SNP_rs.drop_duplicates(subset=['MarkerName'], keep='first')

gwas2018_euro_only = pd.read_csv('/Users/pmerbaum/Desktop/STR_meta/ALS_sumstats_EUR_only.txt.gz', sep = '\t')
only_eur = pd.merge(sorted_SNP_rs_unique, gwas2018_euro_only, left_on='RS', right_on='ID', how='left')
#now only leave genome-wide significant SNPs for a p-p plot
subset_eur = only_eur[((only_eur['P-value_x'] < 5e-08) | (only_eur['P-value_y'] < 5e-08))] 

subset_eur = subset_eur.dropna(subset=['ID'])
subset_eur['-logP_x'] = -np.log10(subset_eur['P-value_x'])
subset_eur['-logP_y'] = -np.log10(subset_eur['P-value_y'])

plt.figure(figsize=(8, 6))
sns.scatterplot(x='-logP_x', y='-logP_y', data=subset_eur)

# Add a trendline
max_value = max(subset_eur['-logP_x'].max(), subset_eur['-logP_y'].max())
plt.plot([0, max_value], [0, max_value], color='gray', linestyle='--')

# Mark significant p-values only in one dataset with red dots
significant_only_in_x = (subset_eur['-logP_x'] < -np.log10(5e-08)) & ~(subset_eur['-logP_y'] < -np.log10(5e-08))
significant_only_in_y = (subset_eur['-logP_y'] < -np.log10(5e-08)) & ~(subset_eur['-logP_x'] < -np.log10(5e-08))

# Plot significant points
plt.scatter(subset_eur.loc[significant_only_in_x, '-logP_x'], 
            subset_eur.loc[significant_only_in_x, '-logP_y'], 
            color='red', label='Significant only in STR GWAS')

plt.scatter(subset_eur.loc[significant_only_in_y, '-logP_x'], 
            subset_eur.loc[significant_only_in_y, '-logP_y'], 
            color='pink', label='Significant only in GWAS 2018')

# Add labels and title
plt.xlabel('-logP STR GWAS')
plt.ylabel('-logP GWAS 2018')
plt.title('P-p plot for STR EHv500_ILL174K_b38_DF2_jan2024 meta (DR>0.7) and GWAS 2018 EUR only')

significant_x_marker_ids = subset_eur.loc[significant_only_in_x, 'Marker_ID']
significant_y_marker_ids = subset_eur.loc[significant_only_in_y, 'Marker_ID']


################################################################################
#comparing STR GWAS and replication study
rep = pd.read_csv('/Users/pmerbaum/Desktop/DF3B_noDF2_nodub_pheno_noNA_panel_EHv500_ILL174K_b38_DF2_jan2024_allCHR_saige_step2_assoc.txt.gz', sep = '\t')
gwas = pd.read_csv('/Users/pmerbaum/Desktop/all_strata_all_chr_genotype_EHv500_ILL174K_b38_DF2_jan2024_meta_STRs.txt.gz', sep = '\t')
rep['threshold'] = rep['Allele1'].str.extract(r'T(\d+)')
gwas['threshold'] = gwas['MarkerName'].str.extract(r'T(\d+)')

merge = pd.merge(gwas, rep, left_on=['CHR', 'POS','threshold'], right_on=['CHR', 'POS','threshold'], how='outer')
significance_threshold = 5e-8

# Filter rows where either P-value or p.value is significant
significant_rows = merge[(merge['P-value'] < significance_threshold) | (merge['p.value'] < significance_threshold)]

# Extract the relevant P-values from the filtered rows (without sorting)
pval1 = significant_rows['P-value']  # P-values from df1 (P-value column)
pval2 = significant_rows['p.value']  # P-values from df2 (p.value column)

# Remove any NaN values
pval1 = pval1.dropna()
pval2 = pval2.dropna()

# Define the conditions for coloring without sorting
significant_only_in_x = (pval1 < significance_threshold) & (pval2 >= significance_threshold)
significant_only_in_y = (pval2 < significance_threshold) & (pval1 >= significance_threshold)

# Plot the PP-plot for significant rows
plt.figure(figsize=(8, 6))

# Plot points where both P-values are significant
plt.scatter(-np.log10(pval1), -np.log10(pval2), color='lightsteelblue', alpha=0.6, label='Both significant')

# Add points where significant only in STR GWAS (significant only in df1)
plt.scatter(-np.log10(pval1[significant_only_in_x]), -np.log10(pval2[significant_only_in_x]), 
            color='red', label='Significant only in STR GWAS')

# Add points where significant only in Replication study (significant only in df2)
plt.scatter(-np.log10(pval1[significant_only_in_y]), -np.log10(pval2[significant_only_in_y]), 
            color='pink', label='Significant only in replication')

# Add a reference line for y = x (perfect match line)
plt.plot([0, max(-np.log10(pval1)), max(-np.log10(pval2))], [0, max(-np.log10(pval1)), max(-np.log10(pval2))], 
         color='gray', linestyle='--', lw=2)  # Gray line

# Labels and title
plt.title('PP-Plot of Significant P-values')
plt.xlabel('-log10(P-value) from STR GWAS')
plt.ylabel('-log10(P-value) from replication study')

# Show the legend
plt.legend()
