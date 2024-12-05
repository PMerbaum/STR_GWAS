#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 09:44:09 2024

@author: pmerbaum
"""
input_file = snakemake.input[0]
output_file = snakemake.output[0]
panel_short = snakemake.params['panel_short']

def process_line(line, chrom_index, pos_index, snpid_index, allele1_index):
    columns = line.strip().split(' ')
    chrom = columns[chrom_index]
    pos = columns[pos_index]
    new_snpid = f'{panel_short}_chr{chrom}_{pos}'

    columns[snpid_index] = new_snpid

    return '  '.join(columns) + '\n'

# Open input file for reading
with open(input_file, 'r') as f_in:
    # Read header
    header = next(f_in).strip()
    header_items = header.split()

    # Find the index of POS and SNPID columns
    pos_index = header_items.index('POS')
    chrom_index = header_items.index('CHR')
    snpid_index = header_items.index('SNPID')
    allele2_index = header_items.index('Allele2')
    allele1_index = header_items.index('Allele1')

    # Open output file for writing
    with open(output_file, 'w') as f_out:
        # Write header to output file
        f_out.write(header + '\n')

        # Process each line in the input file
        for line in f_in:
            updated_line = process_line(line, chrom_index, pos_index, snpid_index, allele1_index)
            # Write updated line to output file
            f_out.write(updated_line)