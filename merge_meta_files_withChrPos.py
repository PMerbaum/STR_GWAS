#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 09:58:47 2024

@author: pmerbaum
"""

import pandas as pd
import os

# Define the list of s and chromosome numbers
s_numbers = [1, 2, 3, 4, 6, 7]
chromosome_numbers = range(1, 23)

# Iterate through each chromosome file
for chr_number in chromosome_numbers:
    # Read the file
    file_path = f"all_strata_panel_EHv500_ILL174K_b38_DF2_jan2024_chr{chr_number}_1_1.txt"
    data = pd.read_csv(file_path, sep='\t')  # Adjust parameters as needed

    # Add CHR column with chromosome number
    data['CHR'] = chr_number

    # Extract base pair position from MarkerName column to create POS column
    data['POS'] = data['MarkerName'].str.extract(r'STR_ILL_174k_(\d+)_T\d+').astype(int)
    data['POS']

    # Save the modified data frame back to a new file
    new_file_path = f"all_strata_panel_EHv500_ILL174K_b38_DF2_jan2024_chr{chr_number}.txt"
    data.to_csv(new_file_path, sep='\t', index=False)
    
    
dfs = []

# Iterate through each chromosome file
for chr_number in range(1, 23):
    # Read the modified file
    file_path = f"all_strata_panel_EHv500_ILL174K_b38_DF2_jan2024_chr{chr_number}.txt"
    data = pd.read_csv(file_path, sep='\t')
    
    # Append the DataFrame to the list
    dfs.append(data)

# Concatenate all DataFrames in the list
combined_df = pd.concat(dfs, ignore_index=True)

# Save the combined DataFrame to a new file
combined_file_path = "all_strata_all_chr_EHv500_ILL174K_b38_DF2_meta.txt"
combined_df.to_csv(combined_file_path, sep='\t', index=False)