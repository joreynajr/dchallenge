#!/usr/bin/env python
# coding: utf-8

# Script adds the distance column to eQTL files that are missing it.

import os
import gzip
import pandas as pd
import re
import argparse
import subprocess

#########################################################################################
# Command line interface
#########################################################################################
parser = argparse.ArgumentParser(description='Add the distance column to the eqtl file.')
parser.add_argument('gencode', help='Path to Gencode file.') 
parser.add_argument('eqtl', help='Path to eQTL file.') 
parser.add_argument('pos', type=int, help='Indicate the column number with the position field.') 
parser.add_argument('geneID', type=int, help='Indicate the column number with the gene id field.') 
parser.add_argument('outfn', help='Path to output eQTL file with a new distance column.') 
parser.add_argument('--header', action='store_true', help='Indicate the presence of a header.') 
params = parser.parse_args()

# load gencode
print("# load gencode")
gencode = pd.read_table(params.gencode, header=None)
gencode.columns = ['chr', 'start', 'end', 'strand', 'type', 'gene_id', 'gene_name']
gencode = gencode.loc[gencode['type'] == 'gene']
print(len(gencode))
gencode = gencode.loc[~gencode.duplicated('gene_id', keep='first')]
print(len(gencode))
gencode.set_index('gene_id', inplace=True)
gencode = gencode[['start', 'end', 'strand']].to_dict('index')
print(len(gencode))


#print(gencode)
print(len(gencode))

eqtl = params.eqtl
outfn = params.outfn

if os.path.splitext(eqtl)[1] == '.gz':
    
    intermfn = params.outfn.replace('.gz', '')

    print('Dealing with .gz files.')
    
    with open(outfn, 'w') as fw:

        # skip the header
        if params.header:
            fw.write('dist\n')

        # calculate distance
        gzip = subprocess.Popen(['gzip', '-cdfq', eqtl], stdout=subprocess.PIPE)
        for i, line in enumerate(gzip.stdout):

            # skip the header
            if i == 0:
                continue

            if i % 1000000 == 0:
                print('\tWorking on line: {}.'.format(i))
            
            # get eqtl info 
            eqtl_info = line.decode().split()

            pos = int(eqtl_info[params.pos - 1])
            gene_info = gencode[eqtl_info[params.geneID - 1]]

            # calculate distance based on the strand information
            if gene_info['strand'] == '+':
                dist = int(abs(gene_info['start'] - pos))
            elif gene_info['strand'] == '-':
                dist = int(abs(gene_info['end'] - pos))
            else:
                raise Exception('Mistake, {} is not a strand orientation.'.format(gene_info.strand))

            fw.write(str(dist) + '\n')
            
else:
    msg = 'Incorrect file types. Both input and output have to use the same extension and '
    msg += 'these extensions must be ".gz". Later developments will included tsv files.'
    raise Exception(msg)
