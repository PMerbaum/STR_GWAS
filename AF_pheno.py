#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 09:45:38 2024

@author: pmerbaum
"""

import pandas as pd
import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#First, extract lengths per sample

dfs = []
directory = '/Users/pmerbaum/Desktop/'

for file_path in glob.glob('STR_length_sig_s*_batch*_chr*.txt'):
    s_number = int(file_path.split('_')[-3][1:])
    batch_number = int(file_path.split('_')[-2][5:])
    chr_number = int(file_path.split('_')[-1].split('.')[0][3:])

    with open(file_path, 'r') as file:
        lines = file.readlines()

        sample_ids = lines[7].strip().split('\t')[9:] 
        str_ids = []
        str_lengths = []
        sample_ids_all = []
        strata = []
        batches = []
        for i in range(8, len(lines), 1):
                    str_id_line = lines[i].strip().split('\t')
                    str_id = str_id_line[2]
                    lengths = str_id_line[9:]
                    for length, sample_id in zip(lengths, sample_ids):
                        str_ids.append(str_id)
                        str_lengths.append(length)
                        sample_ids_all.append(sample_id)
                        strata.append(s_number)
                        batches.append(batch_number)
        df = pd.DataFrame({'STR_ID': str_ids, 'Sample_ID': sample_ids_all, 'Length': str_lengths, 'Strata': strata, 'Batch': batches, 'CHR':chr_number})
        dfs.append(df)
        

sig_length_df = pd.concat(dfs, ignore_index=True)
split_values = sig_length_df['Length'].str.split('/')
swapped_values = split_values.apply(lambda x: f'{x[1]}/{x[0]}' if int(x[0]) > int(x[1]) else '/'.join(x))
sig_length_df['Length'] = swapped_values

split_values = sig_length_df['Length'].str.split('/')
split_values = split_values.apply(lambda x: [int(e) for e in x])
sig_length_df["Short_length"] = split_values.str[0].astype(int)
sig_length_df["Long_length"]= split_values.str[1].astype(int)

###############################################################################
#Second, extract AF per batch

results = []
directory = '/Users/pmerbaum/Desktop/'

# Iterate over all files in the directory
for file_path in glob.glob(directory + 'output_s*_batch*_chr*.txt'):
    # Extract s{number} and batch{number} from the file path
    s_number = int(file_path.split('_')[1][1:])
    batch_number = int(file_path.split('_')[2][5:])
    chr_number = int(file_path.split('_')[-1].split('.')[0][3:])
    with open(file_path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                parts = line.split('\t')
                unit = parts[2].split('_')[-1]
                values = parts[4].split(',')
                amount_units = [v.count(unit) for v in values]
                # Extract values per unit length in eighth column
                values_per_unit = []
                values_per_unit.extend(parts[7].split('AF=')[1].split(','))
                values_per_unit[-1] = values_per_unit[-1][:-4]
                
                # Append the results to the list
                results.append({'STR_ID': parts[2], 'Amount_Units': amount_units, 'Values_Per_Unit': values_per_unit, 'Strata':s_number, 'Batch':batch_number, 'CHR': chr_number})

# Create a DataFrame from the results
result_df = pd.DataFrame(results)
    
##############################################################################

final_df = pd.merge(sig_length_df, result_df, on=['STR_ID', 'Strata', 'Batch'])



##############################################################################

#extract phenotype per sample

strata_numbers = [1,2,3,4,6,7]
directory = '/Users/pmerbaum/Desktop/STR_meta/phenofiles/'

dfs_pheno = []
for pheno in strata_numbers:
    pheno_PC = pd.read_csv(directory + f's{pheno}.no_overlap.hm3_saige_step1.pheno', sep = '\t')
    pheno= pheno_PC[['IID', 'Pheno']].copy()
    pheno.rename(columns={'IID': 'Sample_ID'}, inplace=True)
    pheno['Sample_ID'] = pheno['Sample_ID'].apply(lambda x: x + '_' + x)
    dfs_pheno.append(pheno)
    
phenotype = pd.concat(dfs_pheno, ignore_index=True)

##############################################################################

all_info = pd.merge(final_df, phenotype, on = ['Sample_ID'])

def find_corresponding_value(row):
    long_length = row['Long_length']
    amount_units = row['Amount_Units']
    values_per_unit = row['Values_Per_Unit']
    
    try:
        index = amount_units.index(long_length)
        corresponding_value = values_per_unit[index]
        return corresponding_value
    except ValueError:
        return None


all_info['Corresponding_Value'] = all_info.apply(find_corresponding_value, axis=1)

##############################################################################
#visualize now
all_info = pd.read_csv('/Users/pmerbaum/Desktop/STR_length_AF_pheno_info.txt.gz', sep = '\t')

def AF_length_pheno_str_visualisation(df):
    # Iterate over unique values of 'STR_ID'
    for str_id in df['STR_ID'].unique():
        # Subset the DataFrame for the current 'STR_ID' and drop rows with NA in 'Allele Frequencies'
        subset_df = df[df['STR_ID'] == str_id].dropna(subset=['Corresponding_AF'])
        
        # Create a boxplot for the current subset
        sns.boxplot(x='Long_length', y='Corresponding_AF', hue='Pheno', data=subset_df, showfliers=False)
        plt.title(f'{str_id}: Repeat Lengths vs. Allele Frequencies per phenotype')
        plt.xlabel('Repeat Lengths')
        plt.ylabel('Allele Frequencies')
        plt.legend(title='Pheno', loc = 'upper left')
        plt.show()
        
AF_length_pheno_str_visualisation(all_info)

#plots a frequency table per phenotype for each STR separatelly    
def plot_repeat_length_frequency(df):
    for str_id in df['STR_ID'].unique():
        subset_df = df[df['STR_ID'] == str_id]
        # Create separate DataFrames for each phenotype
        phenotype_0 = subset_df[subset_df['Pheno'] == 0]
        phenotype_1 = subset_df[subset_df['Pheno'] == 1]

        # Calculate normalized counts for both phenotypes
        total_samples = len(subset_df)
        normalized_counts_0 = len(phenotype_0) / total_samples
        normalized_counts_1 = len(phenotype_1) / total_samples

        # Plot histograms with normalized counts for both phenotypes on the same graph
        plt.figure(figsize=(8, 6))

        plt.hist(phenotype_0['Long_length'], bins=10, color='blue', alpha=0.7, label=f'Phenotype 0 (Normalized: {normalized_counts_0:.2f})')
        plt.hist(phenotype_1['Long_length'], bins=10, color='orange', alpha=0.7, label=f'Phenotype 1 (Normalized: {normalized_counts_1:.2f})')

        plt.title(f'Repeat Length Frequency (STR ID: {str_id})')
        plt.xlabel('Repeat Length')
        plt.ylabel('Frequency')
        plt.legend()

        plt.show()

plot_repeat_length_frequency(all_info)

#show distribution of repeat lengths per phenotype for each STR
def repeat_length_percentage_pheno_str_visualisation(df):
    # Iterate over unique values of 'STR_ID'
    for str_id in df['STR_ID'].unique():
        # Subset the DataFrame for the current 'STR_ID'
        subset_df = df[df['STR_ID'] == str_id]

        # Create pivot table to count occurrences of each repeat length for each phenotype
        pivot_table = subset_df.pivot_table(index='Pheno', columns='Short_length', aggfunc='size', fill_value=0)

        # Calculate total number of alleles for each phenotype
        total_alleles_pheno_0 = pivot_table.loc[0].sum()
        total_alleles_pheno_1 = pivot_table.loc[1].sum()

        # Divide count of each repeat length by the total number of alleles to get fraction
        pivot_table.loc[0] /= total_alleles_pheno_0
        pivot_table.loc[1] /= total_alleles_pheno_1

        # Stack the pivot table to have Pheno as a column
        pivot_table = pivot_table.stack().reset_index(name='Fraction')

        # Plot repeat length percentage for each phenotype
        plt.figure(figsize=(10, 5))
        sns.barplot(data=pivot_table, x='Short_length', y='Fraction', hue='Pheno', palette={0: 'lightblue', 1: 'coral'})
        plt.title(f'{str_id}: Repeat Lengths vs. Number of alleles (%)')
        plt.xlabel('Repeat Lengths')
        plt.ylabel('Number of alleles (%)')
        plt.legend(title='Phenotype')
        plt.savefig(f'Short_length_vs_allele_fraction_{str_id}.png' )
        
repeat_length_percentage_pheno_str_visualisation(all_info)



###############################################################################
#extracting imputation information for SNPs per batch
snp_info = []
directory = '/hpc/hers_en/pmerbaum/str_imputation/impute/data/imputed/snps'

for file_path in glob.glob('s*.batch*.no_overlap.hm3_panel_EHv500_ILL174K_b38_DF2_jan2024_chr*_SNP_only_imputed_qual.txt.gz'):
    snp_data = pd.read_csv(file_path, sep ='\t')
    bad_imputed = snp_data.loc[snp_data['DR2'] <= 0.7, 'ID']
    snp_info.append(bad_imputed)
    
snp_info.to_csv('SNPs_dr2_filter_0.7_all_batches.txt.gz', sep = '\t', compression='gzip')
    
    

