#!/usr/bin/env python
# coding: utf-8

# Script adds the distance column to eQTL files that are missing it.

import os
import gzip
import pandas as pd
import re
import argparse

#########################################################################################
# Command line interface
#########################################################################################
parser = argparse.ArgumentParser(description='Add the distance column to the eqtl file.')
parser.add_argument('gencode', help='Path to Gencode file.') 
parser.add_argument('eqtl', help='Path to eQTL file.') 
parser.add_argument('outfn', help='Path to output eQTL file with a new distance column.') 
parser.add_argument('--header', action='store_true', help='Indicate the presence of a header.') 
params = parser.parse_args()

# load gencode
print("# load gencode")
gencode = pd.read_table(params.gencode, header=None)
gencode.columns = ['chr', 'start', 'end', 'strand', 'gene_name', 'score', 'gene_id']
gencode.set_index('gene_id', inplace=True)

eqtl = params.eqtl
outfn = params.outfn

if os.path.splitext(eqtl)[1] == '.gz':
    
    print('Dealing with .gz files.')
    
    with gzip.open(eqtl) as fr, gzip.open(outfn, 'wb') as fw:

        # if the header is present then include it
        if params.header == True:
            # get the header which is composed of the first line
            header = next(fr)

            # getting a list of the current header
            cheader = header.decode().replace('\r', '').strip().split()

            # writing the new header
            nheader = header.decode().replace('\r', '').strip()
            nheader = (nheader + '\tdist\n').encode()
            fw.write(nheader)

        # calculate distance
        for i, line in enumerate(fr):
            
            # get eqtl info 
            eqtl_info = line.decode().replace('\r', '').strip().split()
            eqtl_info = pd.Series(eqtl_info, index=cheader)    
            eqtl_info.position = int(eqtl_info.position)
            gene_info = gencode.loc[eqtl_info.gene_id]

            # calculate distance based on the strand information
            if gene_info.strand == '+':
                dist = abs(gene_info.start - eqtl_info.position)
            elif gene_info.strand == '-':
                dist = abs(gene_info.end - eqtl_info.position)
            else:
                raise Exception('Mistake, {} is not a strand orientation.'.format(gene_info.strand))

            # write out the new eQTL line
            s = eqtl_info.to_string(header=False, index=False)
            s = re.sub('\n\s*', '\t', s)
            s += '\t{}\n'.format(dist)
            fw.write(s.encode())
            
else:
    msg = 'Incorrect file types. Both input and output have to use the same extension and '
    msg += 'these extensions must be ".gz". Later developments will included tsv files.'
    raise Exception(msg)

#elif os.path.splitext(eqtl)[1] == '.tsv':
#    
#    print('Dealing with .tsv files.')
#    
#    # calculate distance
#    with open(eqtl) as fr, open(outfn, 'w') as fw:
#
#        header = next(fr)
#        cheader = header.strip().split()
#        new_header = header + '\tdist\n'
#        fw.write(new_header)
#
#        for line in fr:
#
#            # get eqtl info 
#            eqtl_info = line.strip().split()
#            eqtl_info = pd.Series(eqtl_info, index=cheader)
#            eqtl_info.position = int(eqtl_info.position)
#            gene_info = gencode.loc[eqtl_info.gene_id]
#
#            # calculate distance based on the strand information
#            if gene_info.strand == '+':
#                dist = abs(gene_info.start - eqtl_info.position)
#            elif gene_info.strand == '-':
#                dist = abs(gene_info.end - eqtl_info.position)
#            else:
#                raise Exception('Mistake, {} is not a strand orientation.'.format(gene_info.strand))
#
#            # write out the new eQTL line
#            s = eqtl_info.to_string(header=False, index=False)
#            s = re.sub('\n\s*', '\t', s)
#            s += '\n'
#            fw.write(s)
#            
#else:
#    msg = 'Incorrect file types. Both input and output have to use the same extension and '
#    msg += 'these extensions must either be .gz or .tsv'
#    raise Exception(msg)
