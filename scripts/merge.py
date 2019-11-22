#!/usr/bin/env python
"""Merge STR loci bounds from multiple samples
"""

import argparse
import sys
import os
import collections
import pandas as pd
import pyranges as pr

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Merge overlapping or nearby loci with the same STR repeatunit.')
    parser.add_argument(
        '--bounds', type=str, nargs='+', required = True,
        help='-bounds.txt files from different samples produced by STRling.')
    parser.add_argument(
        '--bed', type=str, nargs='+', required = False,
        help='bed file in the format: chr start end repeatunit [name]')
    parser.add_argument(
        '--out', type=str, default = 'merged.bed',
        help='Output bed file in the format: chr start end repeatunit')
    parser.add_argument(
        '--slop', type=int, default=50,
        help='Merge loci that are within this many bp of each other and have the same repeat unit.')
    return parser.parse_args()

def get_sample(fullpath):
    """Get the sample ID from the filename"""
    basename = os.path.basename(fullpath)
    return(basename.rsplit('-', maxsplit = 1)[0])

def parse_bounds(filename):
    """Parse bounds file produced by STRling"""
    sample_id = get_sample(filename)
    try:
        df = pd.read_csv(filename, sep = '\t', header = None,
                            names = ['Chromosome', 'Start', 'End', 'center_mass', 'n_left', 
                                'n_right', 'n_total', 'repeatunit', 'name', 'median_depth'],
                            usecols = ['Chromosome', 'Start', 'End', 'repeatunit'])
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    df['sample'] = sample_id
    #df = df[['Chromosome', 'Start', 'End', 'repeatunit', 'sample']]
    return(df)

def parse_bed(filename):
    """Parse bed file with the columns: chr start end repeatunit [name]"""
    sample_id = get_sample(filename)
    try:
        df = pd.read_csv(filename, delim_whitespace = True, header = None,
                            names = ['Chromosome', 'Start', 'End', 'repeatunit', 'name'],
                                usecols = ['Chromosome', 'Start', 'End', 'repeatunit'])
    except ValueError:
        df = pd.read_csv(filename, delim_whitespace = True, header = None,
                            names = ['Chromosome', 'Start', 'End', 'repeatunit'])
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    for repeatunit in df['repeatunit']:
        for base in repeatunit:
            if base not in ['A', 'T', 'C', 'G']:
                sys.exit('ERROR: Non-DNA found in the third column of {}: {}\n'.format(filename, repeatunit))
    df['sample'] = sample_id
    return(df)

def merge_subset(genotype_df, repeatunit, slop = 0):
    if genotype_df.shape[0] == 0:
        return(genotype_df)

    if 'locus' in genotype_df.columns:
        genotype_df = genotype_df.drop(['locus', 'Count'], axis = 1)

    genotype_pr = pr.PyRanges(genotype_df)

    merged_pr = genotype_pr.merge(slack = slop, count = True).sort()
    merged_pr.locus = merged_pr.Chromosome.astype(str) + '-' + merged_pr.Start.astype(str) + '-' + merged_pr.End.astype(str) + '-' + repeatunit

    # Annotate genotypes with merged locus
    genotype_pr = genotype_pr.join(merged_pr, suffix = '_x', how = 'left')
    genotype_pr = genotype_pr.drop([x for x in genotype_pr.columns if x.endswith('_x')])

    # Transform back to pandas df and rename position columns back to originals
    genotype_df = genotype_pr.df
    return(genotype_df)

def two_stage_merge(genotype_df, repeatunit, slop = 0):

    # Merge with full slop
    genotype_df = merge_subset(genotype_df, repeatunit, slop = slop)

    # Repeat merge with no slop for duplicated cases
    duplicated_loci_data = genotype_df.loc[genotype_df.duplicated(subset = ['locus', 'sample'], keep=False)]
    # remove duplicated loci from original df
    genotype_df.drop_duplicates(subset = ['locus', 'sample'], keep=False, inplace = True)
    duplicated_loci_merged = merge_subset(duplicated_loci_data, repeatunit, slop = 0)
    genotype_df = genotype_df.append(duplicated_loci_merged, sort = False, ignore_index = True)

    # If duplicate cases remain, give them their original loci
    rows_bool = genotype_df.duplicated(subset = ['locus', 'sample'], keep=False)
    genotype_df.loc[rows_bool, 'locus'] = genotype_df.loc[rows_bool, 'Chromosome'].str.cat(genotype_df.loc[rows_bool, ['Start', 'End', 'repeatunit']].astype(str), sep ="-")

    # Fix loci with no overlapping locus from merging
    rows_bool = genotype_df['Count'] == -1
    genotype_df.loc[rows_bool, 'locus'] = genotype_df.loc[rows_bool, 'Chromosome'].str.cat(genotype_df.loc[rows_bool, ['Start', 'End', 'repeatunit']].astype(str), sep ="-")

    return(genotype_df)

def main():
    # Parse command line arguments
    args = parse_args()

    bounds_files = args.bounds
    bed_files = args.bed
    out_file = args.out
    slop = args.slop

    # Check for duplicate samples
    samples = [get_sample(f) for f in bounds_files]
    duplicate_samples = [item for item, count in collections.Counter(samples).items() if count > 1]
    if len(duplicate_samples) > 0:
        sys.exit('ERROR: duplicate samples found: {}\n'.format(duplicate_samples))

    sys.stderr.write('Processing {} samples\n'.format(len(samples)))

    # Parse data
    bounds = pd.concat( (parse_bounds(f) for f in bounds_files), ignore_index = True)
    if bed_files:
        beds = pd.concat( (parse_bed(f) for f in bed_files), ignore_index = True)
        all_inputs = pd.concat( (bounds, beds), ignore_index = True)
    else:
        all_inputs = bounds

    # Merge
    merged = all_inputs.groupby('repeatunit').apply(lambda x: two_stage_merge(x, repeatunit = x.name, slop = slop))
    unique_loci = set(merged['locus'])

    sys.stderr.write('{} overlapping loci were merged into {}, for a total of {} unique loci after merging\n'.format(merged.loc[merged['Count'] > 1, 'Count'].sum(axis=0),
        merged.loc[merged['Count'] > 1, 'Count'].count(),
        len(unique_loci)))

    # Convert loci back to bed format
    write_data = pd.DataFrame([l.split("-") for l in unique_loci])

    # Write all samples to a single file
    write_data.to_csv(out_file, sep= '\t', index = False, header = False, na_rep='NaN')

if __name__ == '__main__':
    main()
