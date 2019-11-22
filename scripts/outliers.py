#!/usr/bin/env python
"""Read STRling output and look for individuals that are outliers at STR loci
"""

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore", FutureWarning)

    import argparse
    import sys
    import glob
    import os
    #import re
    import numpy as np
    import statsmodels.api as sm
    from scipy.stats import norm
    from statsmodels.sandbox.stats.multicomp import multipletests
    import pandas as pd
    import pyranges as pr

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Read STRling output and look for individuals that are outliers at STR loci')
    parser.add_argument(
        '--genotypes', type=str, nargs='+', required = True,
        help='-genotype.txt files for all samples produced by STRling.')
    parser.add_argument(
        '--unplaced', type=str, nargs='+', required = True,
        help='-unplaced.txt files for all samples produced by STRling. Contains the number of unassigned STR reads for each repeat unit.')
    parser.add_argument(
        '--out', type=str, default = '',
        help='Prefix for all output files (suffix will be STRs.tsv) (default: %(default)s)')
    parser.add_argument(
        '--control', type=str, default='',
        help='Input file for median and standard deviation estimates at each locus from a set of control samples. This file can be produced by this script using the emit option. If this option is not set, all samples in the current batch will be used as controls by default.')
    parser.add_argument(
        '--emit', type=str, default='',
        help='Output file for median and standard deviation estimates at each locus (tsv).')
    parser.add_argument(
        '--slop', type=int, default=50,
        help='Merge loci that are within this many bp of each other and have the same repeat unit.')
    return parser.parse_args()

def get_sample(fullpath):
    """Get the sample ID from the filename"""
    basename = os.path.basename(fullpath)
    return(basename.rsplit('-', maxsplit = 1)[0])

def parse_unplaced(filename):
    """Parse unplaced STR read counts"""
    sample_id = get_sample(filename)
    try:
        unplaced_counts = pd.read_csv(filename, delim_whitespace = True, header = None,
                            names = ['repeatunit', 'unplaced_count'])
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    unplaced_counts['sample'] = sample_id
    unplaced_counts = unplaced_counts[['sample', 'repeatunit', 'unplaced_count']]
    return(unplaced_counts)

def parse_genotypes(filename, min_clips = 5):
    """Parse -genotype.txt file produced by STRling"""
    sample_id = get_sample(filename)
    try:
        genotype_data = pd.read_csv(filename, delim_whitespace = True, header = 0)
        genotype_data.rename(columns={'#chrom': 'Chromosome', 'left': 'Start', 'right': 'End'}, inplace = True)
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    if genotype_data.shape[0] == 0: # Check for file with only header
        sys.exit('ERROR: file {0} contained 0 loci.\n'.format(filename))

    # Filter genotypes
    genotype_data = genotype_data[genotype_data['left_clips'] + genotype_data['right_clips'] >= min_clips]
    genotype_data.reset_index(drop=True, inplace=True)

    # Report number of samples
    sys.stderr.write('Sample: {} Loci: {}\n'.format(sample_id, genotype_data.shape[0]))

    genotype_data['sample'] = sample_id
    return(genotype_data)

def merge_subset(genotype_df, repeatunit, slop = 0):
    if genotype_df.shape[0] == 0:
        return(genotype_df)

    if 'locus' in genotype_df.columns:
        genotype_df = genotype_df.drop(['locus', 'Count'], axis = 1)

    genotype_df = genotype_df.rename(columns={'chrom': 'Chromosome', 'left': 'Start', 'right': 'End'})
    genotype_pr = pr.PyRanges(genotype_df)

    merged_pr = genotype_pr.merge(slack = slop, count = True).sort()
    merged_pr.locus = merged_pr.Chromosome.astype(str) + '-' + merged_pr.Start.astype(str) + '-' + merged_pr.End.astype(str) + '-' + repeatunit

    # Annotate genotypes with merged locus
    genotype_pr = genotype_pr.join(merged_pr, suffix = '_x', how = 'left')
    genotype_pr = genotype_pr.drop([x for x in genotype_pr.columns if x.endswith('_x')])

    # Transform back to pandas df and rename position columns back to originals
    genotype_df = genotype_pr.df
    genotype_df.rename(columns={'Chromosome': 'chrom', 'Start': 'left', 'End': 'right'}, inplace = True)
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
    genotype_df.loc[rows_bool, 'locus'] = genotype_df.loc[rows_bool, 'chrom'].str.cat(genotype_df.loc[rows_bool, ['left', 'right', 'repeatunit']].astype(str), sep ="-")

    # Fix loci with no overlapping locus from merging
    rows_bool = genotype_df['Count'] == -1
    genotype_df.loc[rows_bool, 'locus'] = genotype_df.loc[rows_bool, 'chrom'].str.cat(genotype_df.loc[rows_bool, ['left', 'right', 'repeatunit']].astype(str), sep ="-")

    return(genotype_df)

def parse_controls(control_file):
    """Parse control file with columns locus, median and standard deviation"""

    control_estimates = pd.read_csv(control_file, index_col=0, delim_whitespace = True, header = None)

    # Allow for old style column headings, but change to mu and sd.
    if control_estimates.columns[0] in ['mu', 'median'] and control_estimates.columns[1] in ['sd', 'SD']:
        colnames = list(control_estimates.columns)
        colnames[0:2] = ['mu', 'sd']
        control_estimates.columns = colnames
    else:
        raise ValueError(''.join(["The column names in the control file ",
        "don't look right, expecting columns named median, SD ",
        "or mu, sd. Column names are ", str(list(control_estimates.columns)),
        ". Check the file: ", control_file]))
    return(control_estimates)

def hubers_est(x):
    """Emit Huber's M-estimator median and SD estimates.
    If Huber's fails, emit standard median and NA for sd"""
    huber50 = sm.robust.scale.Huber(maxiter=50)

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)

        try:
            mu, s = huber50(np.array(x))
        except (ValueError, RuntimeWarning):
            mu = np.median(x)
            s = np.nan
    return pd.Series({'mu': mu, 'sd': np.sqrt(s)})

def z_score(x, df):
    """Calculate a z score for each x value, using estimates from a pandas data
    frame with the columns 'mu' and 'sd' and index coressponding to the x values"""
    z = (x.transpose() - df['mu'])/df['sd']
    return z.transpose()

def p_adj_bh(x):
    '''Adjust p values using Benjamini/Hochberg method'''
    return multipletests(x, method='fdr_bh', returnsorted = False)[1]

def main():
    # Parse command line arguments
    args = parse_args()

    base_filename = args.out
    emit_file = args.emit
    control_file = args.control
    genotype_files = args.genotypes
    unplaced_files = args.unplaced
    slop = args.slop
    results_suffix = 'STRs.tsv'

    # Check files exist for all samples
    genotype_ids = set([get_sample(f) for f in genotype_files])
    unplaced_ids = set([get_sample(f) for f in unplaced_files])
    if genotype_ids == unplaced_ids:
        all_samples = genotype_ids
    else:
        all_samples = genotype_ids | unplaced_ids
        missing_samples = (all_samples - genotype_ids) | (all_samples - unplaced_ids)
        sys.exit("ERROR: One or more files are missing for sample(s): " + ' '.join(missing_samples))
    
    sys.stderr.write('Processing {0} samples\n'.format(len(all_samples)))

    if len(all_samples) < 2 and control_file == '':
        sys.stderr.write('WARNING: Only 1 sample and no control file provided, so outlier scores and p-values will not be generated.')

    # Parse unplaced data
    unplaced_data = pd.concat( (parse_unplaced(f) for f in unplaced_files), ignore_index = True)

    unplaced_wide = unplaced_data
    # Fill zeros in unplaced counts
    unplaced_wide = unplaced_data.pivot(index='repeatunit', columns='sample',
                    values='unplaced_count').fillna(0)
    unplaced_wide['repeatunit'] = unplaced_wide.index

    sample_cols = list(set(unplaced_data['sample']))
    unplaced_long = pd.melt(unplaced_wide, id_vars = 'repeatunit',
                            value_vars = sample_cols, value_name = 'unplaced_count',
                            var_name = 'sample')

    # Write unplaced read counts
    unplaced_filename = base_filename + 'unplaced.tsv'
    sys.stderr.write('Writing unplaced counts for all samples to {}\n'.format(unplaced_filename))
    unplaced_long.to_csv(unplaced_filename, sep= '\t', index = False, na_rep='NaN')

    # Parse genotype data
    genotype_data = pd.concat( (parse_genotypes(f) for f in genotype_files), ignore_index = True)

    sys.stderr.write('Merging loci if they are within {} bp of another locus with the same repeat unit\n'.format(slop))
    sys.stderr.write('If merging causes duplicate loci within an individual, merging will be repeated on affected loci without slop.\n')
    sys.stderr.write('Any remaining duplicate loci remaining after all merging will be defined as their original coordinates.\n')

    genotype_data = genotype_data.groupby('repeatunit').apply(lambda x: two_stage_merge(x, repeatunit = x.name, slop = slop))

    sys.stderr.write('{} overlapping loci were merged into {}, for a total of {} unique loci after merging\n'.format(genotype_data.loc[genotype_data['Count'] > 1, 'Count'].sum(axis=0),
        genotype_data.loc[genotype_data['Count'] > 1, 'Count'].count(),
        len(set(genotype_data['locus']))))

#XXX some testing here - remove
    genotype_data.reset_index(drop=True, inplace=True)
#    genotype_data.to_csv("test.csv")
    genotype_data.loc[genotype_data.duplicated(subset = ['locus', 'sample'])].to_csv("test-duplicates.tsv", sep = '\t')

    # Calculate median depth per sample
    sample_depths = genotype_data[['sample','depth']].groupby('sample').median(skipna=True)
    sample_depths['sample'] = sample_depths.index
    sample_depths.to_csv(base_filename + 'depths.tsv', sep= '\t', index = False, na_rep='NaN')

    # Fill zeros in genotype 
    sum_str_wide = genotype_data.pivot(index='locus', columns='sample',
                    values='sum_str_counts').fillna(0)
    sum_str_wide['locus'] = sum_str_wide.index

    sample_cols = list(set(genotype_data['sample']))
    sum_str_long = pd.melt(sum_str_wide, id_vars = 'locus',
                            value_vars = sample_cols, value_name = 'sum_str_counts',
                            var_name = 'sample')

    # Add locus info back in 
    genotype_data = pd.merge(sum_str_long, genotype_data, how='left')
    # Fill zeros in additional columns
    genotype_data[['left', 'right', 'total_reads']] = genotype_data[['left',
                                                    'right', 'total_reads']].fillna(0)

    # Use the median sample depth as a proxy for local depth for loci with zero or NA depth
    genotype_data['depth'] = genotype_data['depth'].replace({0: np.nan})
    genotype_data['depth'] = genotype_data[['sample','depth']].groupby('sample'
                            ).transform(lambda x: x.fillna(x.median(skipna=True)))


    # Normalise STR coverage by median coverage
    factor = 1 # This was 100 in STRetch

    genotype_data['sum_str_log'] = np.log2(factor * (genotype_data['sum_str_counts'] + 1) / genotype_data['depth'])
   
    # For each locus, calculate if that sample is an outlier relative to others
    total_reads_wide = genotype_data.pivot(index='locus', columns='sample', values='sum_str_log')


    sample_depths = genotype_data[['sample', 'depth']].groupby('sample').median()
    # Calculate values for if there were zero reads at a locus in all samples
    null_locus_counts = np.log2(factor * (0 + 1) / sample_depths)
    # Add a null locus that has 0 reads for all individuals (so just uses coverage)
    null_locus_counts_est = hubers_est(null_locus_counts)

    # Calculate a z scores using median and SD estimates from the current set
    # of samples

    # Use Huber's M-estimator to calculate median and SD across all samples
    # for each locus
    locus_estimates = total_reads_wide.apply(hubers_est, axis=1)
    # Where sd is NA, replace with the minimum non-zero sd from all loci
    min_sd = np.min(locus_estimates['sd'][locus_estimates['sd'] > 0])
    locus_estimates['sd'].fillna(min_sd, inplace=True)
    # if sd is 0, replace with min_sd #XXX is this sensible?
    if null_locus_counts_est['sd'] == 0:
        null_locus_counts_est['sd'] = min_sd

    # Save median and SD of all loci to file if requested (for use as a
    # control set for future data sets)
    if emit_file != '':

        locus_estimates.loc['null_locus_counts'] = null_locus_counts_est

        n = len(total_reads_wide.columns)
        locus_estimates['n'] = n

        locus_estimates.to_csv(emit_file, sep= '\t')

    # Calculate z scores using median and SD estimates per locus from a
    # provided control set
    if control_file != '':
        # Parse control file
        control_estimates = parse_controls(control_file)
        # Get a list of all loci in the control file but not the sample data
        control_loci_df = control_estimates.iloc[control_estimates.index != 'null_locus_counts']
        control_loci = [x for x in control_loci_df.index if x not in total_reads_wide.index]

        # Extract and order just those control estimates appearing in the current data
        mu_sd_estimates = control_estimates.reindex(total_reads_wide.index)
        # Fill NaNs with null_locus_counts values
        mu_sd_estimates.fillna(control_estimates.loc['null_locus_counts'],
                                inplace=True)
    else:
        # Extract and order estimates to match the current data
        mu_sd_estimates = locus_estimates.reindex(total_reads_wide.index)

    # calculate z scores
    z = z_score(total_reads_wide, mu_sd_estimates)

    # If a control file is given, effectively add zeros counts at all loci in 
    # controls but not in the samples. 
    # These extra rows will dissapear due to a later merge
    if control_file != '': 
        # Create a total_reads_wide as if all loci have zero counts
        null_total_reads_wide = pd.DataFrame(columns = sample_names, index = control_loci)
        null_total_reads_wide.fillna(null_locus_counts, inplace = True)
        # Caculate z scores
        null_z = z_score(null_total_reads_wide, 
                            control_estimates.reindex(null_total_reads_wide.index))
        loci_with_counts = z.index
        z = z.append(null_z)

    if z.shape[0] == 1:
        ids = z.columns # save index order as data gets sorted
        # Calculate p values based on z scores (one sided)
        z_list = list(z.iloc[0])
        pvals = norm.sf(z_list) # no need to adjust p values if one locus

        # Merge pvals and z scores back into genotype_data
        p_z_df = pd.DataFrame({'sample': ids, 'p_adj': pvals, 'outlier': z_list})
        genotype_data = pd.merge(genotype_data, p_z_df)

    elif z.shape[0] > 1:
        # Calculate p values based on z scores (one sided)
        pvals = z.apply(lambda z_row: [norm.sf(x) for x in z_row], axis=1, 
                    result_type='broadcast') # apply to each row

        if pvals.isnull().values.all(): # Don't bother adjusting p values if all are null
            adj_pvals = pvals
        else:
            # Adjust p values using Benjamini/Hochberg method
            adj_pvals = pvals.apply(p_adj_bh, axis=0) # apply to each column
        
        # Merge pvals and z scores back into genotypes
        adj_pvals['locus'] = adj_pvals.index
        pvals_long = pd.melt(adj_pvals, id_vars = 'locus',
                                value_vars = sample_cols, value_name = 'p_adj',
                                var_name = 'sample')
        genotype_data = pd.merge(genotype_data, pvals_long)
        
        z['locus'] = z.index #important to do this only after p values calculated
        z_long = pd.melt(z, id_vars = 'locus',
                        value_vars = sample_cols, value_name = 'outlier', var_name = 'sample')
        genotype_data = pd.merge(genotype_data, z_long)

    elif z.shape[0] == 0:
        pass #XXX raise error. No input data!

    # Specify output data columns
    write_data = genotype_data[['chrom', 'left', 'right', 'locus',
                                    'sample', 'repeatunit',
                                    'allele1_est', 'allele2_est',
                                    'total_reads', 'spanning_reads', 'spanning_pairs',
                                    'left_clips', 'right_clips', 'unplaced_pairs',
                                    'sum_str_counts', 'sum_str_log', 'depth',
                                    'outlier', 'p_adj',
                                    ]]

    #sort by outlier score then estimated size (bpInsertion), both descending
    write_data = write_data.sort_values(['outlier', 'allele2_est'], ascending=[False, False])
    # Convert outlier and p_adj to numeric type and do some rounding/formatting
    write_data['outlier'] = pd.to_numeric(write_data['outlier'])
    write_data['p_adj'] = [ format(x, '.2g') for x in pd.to_numeric(write_data['p_adj']) ]
    write_data = write_data.round({'outlier': 1, 'sum_str_log': 1,
                                    'sum_str_log': 1})
    int_cols = ['left', 'right', 'total_reads', 'sum_str_counts', 'spanning_reads', 'spanning_pairs',
                'left_clips', 'right_clips', 'unplaced_pairs']
    write_data[int_cols] = write_data[int_cols].astype('Int64')

    # Write individual files for each sample, remove rows where locuscoverage == 0
    samples = set(write_data['sample'])
    for sample in samples:
        sample_filename = base_filename + sample + '.' + results_suffix
        sample_df = write_data.loc[write_data['sample'] == sample]
        sample_df = sample_df.loc[sample_df['total_reads'] != 0.0]
        sample_df.to_csv(sample_filename, sep= '\t', index = False, na_rep='NaN')

    # Write all samples to a single file
    all_filename = base_filename + results_suffix
    write_data.to_csv(all_filename, sep= '\t', index = False, na_rep='NaN')

if __name__ == '__main__':
    main()
