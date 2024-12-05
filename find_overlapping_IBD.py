#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 15:41:38 2024

@author: pmerbaum
Find overlapping IBD segments between samples from the same superpedigree 
"""
import pandas as pd
from intervaltree import IntervalTree

columns = ['SampleID1', 'hapl1', 'SampleID2', 'hapl2', 'CHR','Start_bp','End_bp','LOD']
disease_samples = ['LP6005621-DNA_E07_LP6005621-DNA_E07', 
'LP6005621-DNA_F07_LP6005621-DNA_F07', 
'LP6005621-DNA_B12_LP6005621-DNA_B12', 
'LP6005869-DNA_H07_LP6005869-DNA_H07',
'LP6005621-DNA_H12_LP6005621-DNA_H12',
'LP6008191-DNA_G02_LP6008191-DNA_G02']
control_samples = ['LP6005881-DNA_B04_LP6005881-DNA_B04', 'LP6005879-DNA_F05_LP6005879-DNA_F05']

#for germline output
#path='/Users/pmerbaum/Desktop/grape_germline_IBD.txt'
#columns = ['FID1', 'SampleID1','FID2', 'SampleID2','CHR','Start_bp', 'End_bp', 'Start_snp', 'End_snp', 'total_snp','length','Units','Mismatch','S1_hom','S2_hom']

disease_samples = ['0_LP6005621-DNA_E07_LP6005621-DNA_E07', 
'0_LP6005621-DNA_F07_LP6005621-DNA_F07', 
'0_LP6005621-DNA_B12_LP6005621-DNA_B12', 
'0_LP6005869-DNA_H07_LP6005869-DNA_H07',
'0_LP6005621-DNA_H12_LP6005621-DNA_H12',
'0_LP6008191-DNA_G02_LP6008191-DNA_G02']
control_samples = ['0_LP6005881-DNA_B04_LP6005881-DNA_B04', '0_LP6005879-DNA_F05_LP6005879-DNA_F05']

# Function to determine the status of the pair
def determine_status(sample1, sample2):
    if sample1 in disease_samples and sample2 in disease_samples:
        return 'disease-disease'
    elif sample1 in control_samples and sample2 in control_samples:
        return 'control-control'
    elif sample1 in disease_samples and sample2 in control_samples:
        return 'disease-control'
    else:
        return 'unknown'
    
def compute_interval_length(start, end):
    return end - start
    
def create_interval_trees(df):
    trees = {}
    for chrom, group in df.groupby('CHR'):
        tree = IntervalTree()
        for _, row in group.iterrows():
            tree[row['Start_bp']:row['End_bp']] = (row['SampleID1'], row['SampleID2'])
        trees[chrom] = tree
    return trees

def find_overlapping_segments_with_new_positions(tree, other_trees):
    segments = []
    for interval in tree:
        # Check for overlaps with the disease-disease tree
        overlaps = tree.overlap(interval.begin, interval.end)
        if not overlaps:
            continue  # Skip if there are no overlaps

        overlap_start = max([o.begin for o in overlaps])
        overlap_end = min([o.end for o in overlaps])

        # Ensure a valid overlap
        if overlap_start < overlap_end:
            # Check if this overlap exists in other trees
            is_excluded = False
            for other_tree in other_trees:
                other_overlaps = other_tree.overlap(overlap_start, overlap_end)
                if other_overlaps:
                    is_excluded = True
                    break  # Exit loop if any overlap is found in other trees
            
            if not is_excluded:
                length = compute_interval_length(overlap_start, overlap_end)
                segments.append((interval.data[0], interval.data[1], overlap_start, overlap_end, length))
    return segments

######Modifiable part#######
path = '/Users/pmerbaum/Desktop/DF2_NL_b37_chr7.ibd.gz'
pairs = pd.read_csv(path, sep=r'[ \t]+', header=None, names=columns)
pairs[['SampleID1', 'SampleID2']] = pairs[['SampleID1', 'SampleID2']].apply(sorted, axis=1, result_type='expand')


# Apply the function to each row to create the 'status' column
pairs['Status'] = pairs.apply(lambda row: determine_status(row['SampleID1'], row['SampleID2']), axis=1)
pairs['Status'].value_counts()

disease_disease_pairs = pairs[pairs['Status'] == 'disease-disease']
disease_control_pairs = pairs[pairs['Status'] == 'disease-control']
unknown_pais = pairs[pairs['Status'] == 'unknown']
control_control_pairs = pairs[pairs['Status'] == 'control-control']

#change to the df of interest
interval_trees_disease = create_interval_trees(disease_disease_pairs)
interval_trees_dc = create_interval_trees(disease_control_pairs)
interval_trees_cc = create_interval_trees(control_control_pairs)
interval_trees_unknown = create_interval_trees(unknown_pais)

overlapping_segments = []
other_trees =  list(interval_trees_cc.values()) + list(interval_trees_dc.values()) #+ list(interval_trees_unknown.values())

for chrom, tree in interval_trees_disease.items():
    segments = find_overlapping_segments_with_new_positions(tree, other_trees)
    overlapping_segments.extend([(chrom, *segment) for segment in segments])
      

overlap_df = pd.DataFrame(overlapping_segments, columns=['CHR', 'SampleID1', 'SampleID2', 'Start_bp', 'End_bp', 'Length'])
overlap_df = overlap_df.drop_duplicates()
overlap_df = overlap_df.sort_values(by='Length', ascending=False)

overlap_counts = overlap_df.groupby(['CHR', 'Start_bp', 'End_bp']).size().reset_index(name='Count')
max_count = overlap_counts['Count'].max()
print(max_count)

# Filter segments that have the maximum count of overlaps
max_overlapping_segments = overlap_counts[overlap_counts['Count'] == max_count]
print(max_overlapping_segments)

# Extract SampleID1 and SampleID2 columns
max_overlapping_pairs = pd.merge(max_overlapping_segments, overlap_df, on=['CHR', 'Start_bp', 'End_bp'])
pairs_with_max_overlap = max_overlapping_pairs[['SampleID1', 'SampleID2']]
print(set(pairs_with_max_overlap['SampleID1'])) #no 0_LP6008191-DNA_G02_LP6008191-DNA_G02

start_end_tuples = list(zip(overlap_df['Start_bp'], overlap_df['End_bp']))
unique_start_end_tuples = list(set(start_end_tuples))
# Sort the list of tuples by the 'Start_bp' (first element in the tuple)
sorted_start_end_tuples = sorted(unique_start_end_tuples, key=lambda x: x[0])
