#!/usr/bin/env python
# coding: utf-8

import os
import gzip
import argparse 
from statsmodels.sandbox.stats.multicomp import multipletests

#########################################################################################
# Command line interface
#########################################################################################
parser = argparse.ArgumentParser(description="Calculate the fdr")
parser.add_argument('eqtl', help='Path to the eQTL table.')
parser.add_argument('outfn', help='Path to the eQTL table.')
parser.add_argument('pval_col', metavar='pval-col', type=int, help='1-based column index of the p-value.')
parser.add_argument('--in-header',
                    action='store_true',
                    help='Indicate whether we should parse a header in the eqtl file.')
parser.add_argument('--out-header',
                    action='store_true',
                    help='Indicate whether we should parse a output a header in the out file.')
params = parser.parse_args()

print('The command line interface captured the following values.')
print(params)

#########################################################################################
# Get the pvalues and calculated fdr values
#########################################################################################
pvals = []

# find which open function to use 
print("# find which open function to use")
if os.path.splitext(params.eqtl)[1] == '.gz':
    fopen = lambda x: gzip.open(x, 'rt')
else:
    fopen = lambda x: open(x, 'r')

# parse the eqtl file and extact the p-values
print("# parse the eqtl file and extact the p-values")
with fopen(params.eqtl) as fr:

    # parse the header
    if params.in_header:
        next(fr)

    # parse the rest of the file
    for i, line in enumerate(fr):
        
        if (i % 1000000 == 0):
            print('parsing line: {}'.format(i))
            
        pval = float(line.split()[params.pval_col - 1])
        pvals.append(pval)

# Calculate the FDR Benjamini-Hochberg
print("# Calculate the FDR Benjamini-Hochberg")
passed, fdr_vals, alphacSidak, alphacBonf = multipletests(pvals, alpha=0.05, method='fdr_bh')

#########################################################################################
# Get the fdr values as a single column txtg
#########################################################################################

# Save the fdr txt file
print("# Save the fdr txt file")
with open(params.outfn, 'w') as fw:

    # parse the header
    if params.out_header:
        fw.write('fdr\n')

    # write the rest of the file
    for v in fdr_vals:
        fw.write('{}\n'.format(v))

