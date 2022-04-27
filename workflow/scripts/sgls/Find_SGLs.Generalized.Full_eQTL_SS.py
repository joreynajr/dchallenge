#!/usr/bin/env python
# coding: utf-8
import os 
import sys
import pybedtools as pbt
import pandas as pd
import numpy as np
import subprocess as sp
import json
import argparse
from itertools import chain
pbt.helpers.set_tempdir('/mnt/BioHome/jreyna/tmp/')
pbt.set_bedtools_path('/mnt/BioApps/bedtools/bin/')
bgzip = '/mnt/BioApps/tabix/tabix-0.2.6/bgzip'
tabix = '/mnt/BioApps/tabix/tabix-0.2.6/tabix'

#pd.options.display.max_columns = 1000
#pd.options.display.max_rows = 1000

#os.chdir('/mnt/BioHome/jreyna/jreyna/projects/dchallenge/')

# setting basical column names for bedpe
bedpe_6cols = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB']
bedpe_10cols = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB', 'name', 'score', 'strandA', 'strandB']

# # Making a command line interface
parser = argparse.ArgumentParser()
parser.add_argument('--eqtls', type=str, required=True)
parser.add_argument('--colocs', type=str, required=True)
parser.add_argument('--loops', type=str, required=True)
parser.add_argument('--gs', type=str, required=True)
parser.add_argument('--gene-ref', type=str, required=True)
parser.add_argument('--loop-slop', type=int, default=5000, required=False)
parser.add_argument('--outdir', type=str, required=True)

# argumnts to extract the correct fields from the eqtl table 
parser.add_argument('--eqtl-chr', type=int, required=True)
parser.add_argument('--eqtl-pos', type=int, required=True)
parser.add_argument('--eqtl-pvalue', type=int, required=True)
parser.add_argument('--eqtl-beta', type=int, required=True)
parser.add_argument('--eqtl-geneid', type=int, required=True)
parser.add_argument('--eqtl-dist', type=int, required=True)
parser.add_argument('--eqtl-fdr', type=int, required=True)
parser.add_argument('--eqtl-prefix', type=str, default='eqtl')

# argumnts to extract the correct fields from the coloc table
parser.add_argument('--coloc-chr', type=int, default=1)
parser.add_argument('--coloc-pos', type=int, default=2)
parser.add_argument('--coloc-gid', type=int, required=True)
parser.add_argument('--coloc-prefix', type=str, default='coloc')

# argumnts to extract the correct fields from the coloc table
parser.add_argument('--loop-chrA', type=int, default=1)
parser.add_argument('--loop-startA', type=int, default=2)
parser.add_argument('--loop-endA', type=int, default=3)
parser.add_argument('--loop-chrB', type=int, default=4)
parser.add_argument('--loop-startB', type=int, default=5)
parser.add_argument('--loop-endB', type=int, default=6)
parser.add_argument('--loop-score', type=int, default=7);
parser.add_argument('--loop-prefix', type=str, default='loop')

# creating dummy values for debugging, necessary since I started this 
# code within a jupyter notebook
debug = False
if debug == True: 
    debug = []
    debug.append('--eqtls A')
    debug.append('--colocs B')
    debug.append('--loops C')
    debug.append('--gs D')
    debug.append('--gene-ref E')
    debug.append('--loop-slop 5000')
    debug.append('--outdir F')

    debug.append('--eqtl-chr 1')
    debug.append('--eqtl-pos 2')
    debug.append('--eqtl-pvalue 6')
    debug.append('--eqtl-beta 5')
    debug.append('--eqtl-geneid 3')
    debug.append('--eqtl-dist 4')
    debug.append('--eqtl-fdr 7')

    debug.append('--coloc-chr 1')
    debug.append('--coloc-pos 2')
    debug.append('--coloc-gid 10')

    debug.append('--loop-chrA 1')
    debug.append('--loop-startA 2')
    debug.append('--loop-endA 3')
    debug.append('--loop-chrB 4')
    debug.append('--loop-startB 5')
    debug.append('--loop-endB 6')
    debug.append('--loop-score 20')

    # split and chain everything together
    debug = [x.split() for x in debug]
    debug = list(chain(*debug))

    params = parser.parse_args(debug)
    
    params.eqtls = 'results/main/GRCh37/sgls/ImmuNexUT/T1D_34012112_Gaulton/ImmuNexUT/DN_B/eqtls.coloc_filtered.tsv.gz'
    params.colocs = 'results/main/GRCh37/coloc/eQTL_Catalogue/T1D_34012112_Gaulton/ImmuNexUT/DN_B/ldpairs/coloc_ld_snps.txt'
    params.loops = 'results/main/h3k27ac_hichip/CD8_T-cell_naive/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed'
    params.gs = 'results/refs/hg19/hg19.chrom.sizes'
    params.gene_ref = 'results/refs/gencode/v30/gencode.v30.annotation.grch37.bed'
    params.outdir = 'results/main/GRCh37/sgls/ldpairs/ImmuNexUT/T1D_34012112_Gaulton/ImmuNexUT/DN_B/monocyte_naive/script_version/'
else:
    params = parser.parse_args()

# - [X] CD16p_Mono
# - [X] CL_Mono
# - [X] CM_CD8 (T-cell based)
# - [X] DN_B (T-cell based; Good number of overlaps!)
# - [ ] Next
# - [ ] Next
# - [ ] Next

# # Updating column index fields
def get_compliment_cols(curr_col_idxs, full_idxs):
    """
    Function to help replicate R's -c(1,2,3) indexing trick.
    """
    
    full = set(full_idxs)
    curr = set(curr_col_idxs)
    return(sorted(full.difference(curr)))

# convert from 1-based to 0-based indexing and
# store the indexes in lists
eqtl_cols = ['eqtl-chr', 'eqtl-pos', 'eqtl-geneid', 'eqtl-beta', 'eqtl-pvalue', 'eqtl-fdr', 'eqtl-dist']
coloc_cols = ['coloc-chr', 'coloc-pos', 'coloc-gid']
loop_cols = ['loop-chrA', 'loop-startA', 'loop-endA', 'loop-chrB', 'loop-startB', 'loop-endB', 'loop-score'] 

eqtl_cols_indices = []
for k in eqtl_cols:
    params.__dict__[k.replace('-', '_')] -= 1
    eqtl_cols_indices.append(params.__dict__[k.replace('-', '_')])

coloc_cols_indices = []
for k in coloc_cols:
    params.__dict__[k.replace('-', '_')] -= 1
    coloc_cols_indices.append(params.__dict__[k.replace('-', '_')])

loop_cols_indices = []
for k in loop_cols:
    params.__dict__[k.replace('-', '_')] -= 1
    loop_cols_indices.append(params.__dict__[k.replace('-', '_')])

# # Load the colocalization data
coloc = pd.read_table(params.colocs)

# ensuring geneid column is properly names
coloc.rename(columns={coloc.columns[params.coloc_gid]: 'geneid'}, inplace=True)

# make sure gene IDs are not versioned
coloc.iloc[:, params.coloc_gid] = coloc.iloc[:, params.coloc_gid].str.replace('.[0-9]*$', '', regex=True)

# add and sid 
if 'sid' not in coloc.columns:
    coloc.loc[:, 'sid'] = coloc.loc[:, 'chr'].str.replace('chr', '', regex=True)
    coloc.loc[:, 'sid'] += ':' + coloc.loc[:, 'pos'].astype(str)

# reorder the columns so the main coloc columns are in front,
# while the other columns are in back
other_coloc_cols_indices = get_compliment_cols(coloc_cols_indices, np.arange(coloc.shape[1]))
coloc = coloc.iloc[:, coloc_cols_indices + other_coloc_cols_indices]

# remove duplicates snp-gene pairs
# it won't happen often but sometimes the loci we
# choose can overlap and a SNP-gene pair will be retested 
coloc_sig_df = coloc.loc[~coloc.duplicated(subset=['sid', 'geneid'])]

# keeping a copy of this coloc_sig_df
coloc_sig_full = coloc_sig_df.copy(deep=True)

# make a coloc pybed tools object
coloc_sig_df.rename(columns={'pos': 'end'}, inplace=True)
coloc_sig_df.loc[:, 'start'] = coloc_sig_df.loc[:, 'end'] - 1
coloc_sig_df = coloc_sig_df[['chr', 'start', 'end', 'sid']]
coloc_sig_pbt = pbt.BedTool.from_dataframe(coloc_sig_df.iloc[:, 0:4]).sort()

print('There are {} colocalized SNP-gene pairs'.format(coloc_sig_df.shape[0]))

# # Load the gencode coords
# Load the gene data
print('# Load the gene data')

cols = ['chrom', 'start', 'end', 'strand', 'type', 'geneid', 'genename']
gencode = pd.read_table(params.gene_ref, header=None, names=cols)

# extract just the genes
genes_df = gencode.loc[gencode.type.isin(['gene'])]
genes_df = genes_df.loc[~genes_df.duplicated(subset='geneid'), :]
genes_df.loc[:, 'chrom'] = genes_df['chrom'].astype(str)
genes_df = genes_df.iloc[:, [0,1,2,6,5,3]]

# create a copy of the original gene bed before coordinate shrinking
orig_genes_df = genes_df.copy()

# convert the start/end position into start/end for the TSS
# if the gene is + then the start is uses as the tss otherwise
# the end is used as the tss
genes_df.loc[(genes_df.strand == '+'), 'end'] = genes_df.loc[(genes_df.strand == '+'), 'start']
genes_df.loc[(genes_df.strand == '+'), 'start'] = genes_df.loc[(genes_df.strand == '+'), 'start'] - 1
genes_df.loc[(genes_df.strand == '-'), 'end'] = genes_df.loc[(genes_df.strand == '-'), 'end']
genes_df.loc[(genes_df.strand == '-'), 'start'] = genes_df.loc[(genes_df.strand == '-'), 'end'] - 1

# make a genes pbt for intersection
print("# make a genes pbt for intersection")
print(genes_df.head())
genes_pbt = pbt.BedTool.from_dataframe(genes_df).sort()

print('\nThere are {} genes in this GTF-derived file.'.format(genes_df.shape[0]))

# ## Find all genes +/- 1mb of a colocalized SNP
print('# Find all genes +/- 1mb of a colocalized SNP')
slop = 1e6

# get a list of gene names within +- 1mb of the SNPs
coloc_snp_windows = coloc_sig_pbt.slop(b=slop, g=params.gs)
cd_sgls = coloc_snp_windows.intersect(genes_pbt, wa=True, wb=True).to_dataframe()

# reorder for downstream loop intersections 
cd_sgls = cd_sgls.iloc[:, [0,1,2,4,5,6,3,7,8,9]]
cd_sgls.columns = bedpe_6cols + ['sid', 'genename', 'geneid', 'strandB']

# add remaining bedpe-10 columns and reordering
cd_sgls['strandA'] = '+'
cd_sgls['name'] = cd_sgls['sid'] + '-' + cd_sgls['geneid']
cd_sgls['score'] = '.'
new_order = bedpe_10cols + ['sid', 'genename', 'geneid']
cd_sgls = cd_sgls[new_order]

# remove the sloop added in window step
# also ensuring that integer columns remain so, some 
# calculations with these columns can change them into floats
# and cause downstream effects like missing loop intersections
# which found out about 2022-04-23.
cd_sgls['startA'] = (cd_sgls['startA'] + slop).astype(int)
cd_sgls['endA'] = (cd_sgls['endA'] - slop).astype(int)
cd_sgls['startB'] = cd_sgls['startA'].astype(int)
cd_sgls['endB'] = cd_sgls['endB'].astype(int)

# print a quick coloc summary
print('There are {} colocalized snp-gene pairs within +/- 1mb.'.format(cd_sgls.shape[0]))

# # Find the closest gene
print('# Find the closest gene')
closest_gene = coloc_sig_pbt.closest(genes_pbt, d=True, t='all').to_dataframe()
closest_gene = closest_gene.iloc[:, [0,1,2,4,5,6,3,7,8,9,10]]
closest_gene.columns = bedpe_6cols + ['sid', 'genename', 'geneid', 'strand', 'dist']
closest_gene.set_index(['sid', 'geneid'], inplace=True)

# # Load the loop data
print('# Load the H3K27ac HiChIP loops')

print('# Get the loops')
# load the loop data
full_loops = pd.read_table(params.loops)

# reorder the columns so the main loops columns are in front,
# while the other columns are in back
other_loop_cols_indices = get_compliment_cols(loop_cols_indices, np.arange(full_loops.shape[1]))
full_loops = full_loops.iloc[:, loop_cols_indices + other_loop_cols_indices]
full_loops.columns = bedpe_6cols + ['score'] + full_loops.columns.tolist()[7:]

# renaming the loop columns with a loop prefix
rename_cols = ['{}.{}'.format(params.loop_prefix, x) for x in full_loops.columns.tolist()]
full_loops.columns = rename_cols

# add a loop ID for downstream matching
loop_ids = []
for vec in full_loops.values:
    fields = [params.loop_chrA, params.loop_startA, params.loop_endA, 
              params.loop_chrB, params.loop_startB, params.loop_endB]
    loop_id = '{}.{}.{}.{}.{}.{}'.format(*vec[fields])
    loop_ids.append(loop_id)
full_loops.loc[:, 'loop.id'] = loop_ids

# creating a temp loop obj for looping coordinates + loop ID
tmp_loops = full_loops.iloc[:, 0:6]

# applying slop to loops
tmp_loops.iloc[:, params.loop_startA] -= params.loop_slop
tmp_loops.iloc[:, params.loop_endA] += params.loop_slop
tmp_loops.iloc[:, params.loop_startB] -= params.loop_slop
tmp_loops.iloc[:, params.loop_endB] += params.loop_slop

# adding other bed columns
tmp_loops.loc[:, 'loop.name'] = full_loops.loc[:, 'loop.id']
tmp_loops.loc[:, 'loop.score'] = full_loops.iloc[:, 6]
tmp_loops.loc[:, 'loop.strand1'] = '.'
tmp_loops.loc[:, 'loop.strand2'] = '.'

# create a bedtool object
loops = pbt.BedTool.from_dataframe(tmp_loops)
print(tmp_loops.head())
print()
print('There are {} loops in total.'.format(tmp_loops.shape[0]))


# ## Find SNP-Gene overlapping a loop (both anchors)
print('### Find SNP-Gene overlapping a loop (both anchors)')

# re-arranging to fit bedpe format
pr_sgls = cd_sgls.copy()

# loading into pbt
pr_sgls = pbt.BedTool.from_dataframe(pr_sgls)
pr_sgls = pr_sgls.pair_to_pair(loops, type='both', **{'is':True})
pr_sgls = pr_sgls.to_dataframe(disable_auto_names=True, header=None)

if len(pr_sgls) > 0:
    pr_sgls_set = pr_sgls.iloc[:, [10,12]]
    pr_sgls_uniq = set([tuple(x) for x in pr_sgls_set.values.tolist()])
else:
    print('WARNING: Found no overlap between the genes and loops.')
    pr_sgls_uniq = set()

print('There are {} SNP-Gene pairs with a loop.'.format(len(pr_sgls_uniq)))

# # Extract coloc-SNP overlapping an anchor
print('### Find coloc-SNP overlapping an anchor')
if len(pr_sgls) > 0:
    colocSNP_anchors_set = set(pr_sgls.iloc[:, 10].values)
else:
    print('WARNING: Found no overlap between coloc SNPs and loop anchors.')
    colocSNP_anchors_set = set()

print('There are {} SNPs which overlap a loop anchor.'.format(len(colocSNP_anchors_set)))

# # Find coloc-gene overlapping an anchor
print('### Find coloc-gene overlapping an anchor')

if len(pr_sgls) > 0:
    colocGene_anchors_set = pr_sgls.iloc[:, 12].unique()
else:
    print('WARNING: Found no overlap between coloc Genes and loop anchors.')
    colocGene_anchors_set = set()

print('There are {} Genes which overlap a loop anchor.'.format(len(colocGene_anchors_set)))

# ### Creating a mapper between coloc id and loop id
coloc_loop_mapper = pr_sgls.iloc[:, [6, 19]]
coloc_loop_mapper.columns = ['coloc.id', 'loop.id']
coloc_loop_mapper.drop_duplicates(inplace=True)

# # Construct master table
print('# Construct master table')

# begin making the master
master = cd_sgls.copy()

print('Master is starting with {} candidate snp-gene pairs.'.format(master.shape[0]))


# ## Add eQTL fields/details 
print('#### Add eqtl results')

# get SNP eQTL's; filtering required for eQTL Catalogue data
eqtls = pd.read_table(params.eqtls, header=0)

# gathering a dictionary of eqtl renames
# renaming must be one AFTER reordering columns
eqtls_rename_cols = {
    eqtls.columns[params.eqtl_chr]: 'eqtl.chr',
    eqtls.columns[params.eqtl_pos]: 'eqtl.pos',
    eqtls.columns[params.eqtl_geneid]: 'geneid',
    eqtls.columns[params.eqtl_beta]: 'eqtl.beta',
    eqtls.columns[params.eqtl_pvalue]: 'eqtl.pval',
    eqtls.columns[params.eqtl_dist]: 'eqtl.dist', 
    eqtls.columns[params.eqtl_fdr]: 'eqtl.fdr'}

# reorder the columns so the main coloc columns are in front,
# while the other columns are in back. This reordering MUST
# come before the reordering of columns
other_eqtl_cols_indices = get_compliment_cols(eqtl_cols_indices, np.arange(eqtls.shape[1]))

# adding eqtl prefix to the other columns
for i in other_eqtl_cols_indices:
    eqtls_rename_cols.update({eqtls.columns[i]: '{}.{}'.format(params.eqtl_prefix, eqtls.columns[i])})

eqtls = eqtls.iloc[:, eqtl_cols_indices + other_eqtl_cols_indices]

# renaming the main eqtl columns
eqtls.rename(columns=eqtls_rename_cols, inplace=True)

# remove eqtls where chr and pos are empty
# this came up with ImmuNexUT datasets
eqtls = eqtls.loc[~eqtls['eqtl.pos'].isna()]

# make sure the pos column is designated as int
eqtls.loc[:, 'eqtl.pos'] = eqtls['eqtl.pos'].astype(int)

# remove versioned portion of gene ID
eqtls.geneid = eqtls.geneid.str.replace('\.[0-9]*$', '', regex=True)

print(eqtls.head())

eqtls['sid'] = eqtls['eqtl.chr'].replace('chr', '', regex=True).astype(str) + ':' + eqtls['eqtl.pos'].astype(str)
eqtls['flt.is_eqtl_pair'] = 1 # add column to filter on eqtl snp status
print('There are {} eQTLs.'.format(eqtls.shape[0]))

# need to use outer or else you exclude some eQTL's
master = master.merge(eqtls, on=['sid', 'geneid'], how='left')

# update eqtl snp status
master.loc[master['flt.is_eqtl_pair'].isna(), 'flt.is_eqtl_pair' ] = 0 
master.loc[:, 'flt.is_eqtl_pair'] = master.loc[:, 'flt.is_eqtl_pair'].astype(int)

# add gene names to entries with a missing name (after adding eQTL info)
master.loc[master.genename.isna(), 'genename'] = master.loc[master.genename.isna(), 'geneid']

# add missing chrA, chrB, startA and endA data for the eQTL rows
master.loc[master.chrA.isna(), 'chrA'] = 'chr' + master.loc[master.chrA.isna(), 'sid'].str.replace(':[0-9]+', '', regex=True)
master.loc[master.chrB.isna(), 'chrB'] = 'chr' + master.loc[master.chrB.isna(), 'sid'].str.replace(':[0-9]+', '', regex=True)
master.loc[master.startA.isna(), 'startA'] = (master.loc[master.startA.isna(), 'sid'].str.replace('[0-9]+:', '', regex=True)).astype(int)
master.loc[master.startA.isna(), 'startA'] -= 1     
master.loc[master.endA.isna(), 'endA'] = master.loc[master.endA.isna(), 'sid'].str.replace('[0-9]+:', '', regex=True)

# add indicator of eqtl gene presence, it is often the case that a gene is an egene 
# but it's eSNP is not a colocSNP, so it's helpful to have this column
master.loc[:, 'flt.is_eqtl_gene'] = (master.geneid.isin(eqtls.geneid.unique())).astype(int)

print('After outer merging with eqtls, master has {} snp-gene pairs.'.format(master.shape[0]))


# ## Add gene meta fields/columns 
print('### Add gene meta data')

# add back the original gene start and end
master = master.merge(orig_genes_df[['start', 'end', 'geneid', 'strand']],
                      on='geneid', how='left')

# convert the startB/endB position into startB/endB for the TSS
# if the gene is + then the startB is uses as the tss otherwise
# the endB is used as the tss
master.loc[(master.strand == '+'), 'endB'] = master.loc[(master.strand == '+'), 'start']
master.loc[(master.strand == '+'), 'startB'] = master.loc[(master.strand == '+'), 'start'] - 1

master.loc[(master.strand == '-'), 'endB'] = master.loc[(master.strand == '-'), 'end']
master.loc[(master.strand == '-'), 'startB'] = master.loc[(master.strand == '-'), 'end'] - 1

# convert the coordinates from floats to ints
master.startA = master.startA.astype(int)
master.startB = master.startB.astype(int)
master.endA = master.endA.astype(int)
master.endB = master.endB.astype(int)
master.rename(columns={'start': 'gene.start', 'end': 'gene.end', 'strand': 'gene.strand'}, inplace=True)

# # Add boolean filter about closests gene
print('# # Add boolean filter about closests gene')

# check for the closets gene
closets_check = [0] * master.shape[0]
for i, sr in master.iterrows():
    rs_gene = (sr.sid, sr.geneid)
    if rs_gene in closest_gene.index:
        closets_check[i] = 1
master['flt.is_closest_gene'] = closets_check
print('Found {} closest genes.'.format(sum(master['flt.is_closest_gene'])))


# ## Add colocalization data
print('### Add colocalization data')

# rename the columns with coloc prefix
rename_cols = []
for col in coloc_sig_full:
    if col not in ['sid', 'geneid']:
        rename_cols.append('{}.{}'.format(params.coloc_prefix, col))
    else:
        rename_cols.append(col)

# add colocalization data for SNP and is_coloc_snp columns
tmp_coloc = coloc_sig_full.copy()
tmp_coloc.columns = rename_cols
master = master.merge(tmp_coloc, on=['sid', 'geneid'], how='left')

# ### Add column to filter on coloc snp status

print('### Add column to filter on coloc snp status')
master['flt.is_coloc_pair'] = (~master['coloc.pp_H4_Coloc_Summary'].isna()).astype(int)

print('After left merging master with the colocalization table there are {} entries.'.format(master.shape[0]))
after_cnt = sum(master['flt.is_coloc_pair'])
print('Checking if I have the correct number of colocalized SNPs: before {}; after {}.'.      format(coloc_sig_df.shape[0], after_cnt))

# ## Add loop data

# add loop meta-data
master = master.merge(coloc_loop_mapper, left_on='name', right_on='coloc.id', how='left')
master = master.merge(full_loops, on='loop.id', how='left')

# add boolean field for has_fithichip_loop
# check for the loop gene
loop_check = [0] * master.shape[0]
for i, sr in master.iterrows():
    # check closest gene
    rs_gene = (sr.sid, sr.geneid)
    if rs_gene in pr_sgls_uniq:
        loop_check[i] = 1       

master['flt.has_fithichip_loop'] = loop_check
print('There are {} SNP-Gene loops.'.format(sum(loop_check)))
print()

# add boolean field for has_snp_anchor
print('Add coloc-snp boolean column')
master['flt.has_snp_anchor'] = master.sid.isin(colocSNP_anchors_set).astype(int)
print('There are {} coloc-SNP anchors.'.format(sum(master['flt.has_snp_anchor'])))

master['flt.has_gene_anchor'] = master.geneid.isin(colocGene_anchors_set).astype(int)
print('There are {} coloc-Gene anchors.'.format(sum(master['flt.has_gene_anchor'])))

# adding a convenience columns for SGL analysis with 
# SNP-gene pairs coming from colocalization and eqtl's
master['flt.is_coloc_sgl'] = ((master['flt.has_fithichip_loop'] == 1) & (master['flt.is_coloc_pair'] == 1)).astype(int)
master['flt.is_eqtl_sgl'] = ((master['flt.has_fithichip_loop'] == 1) & (master['flt.is_eqtl_pair'] == 1)).astype(int)

# ## Do the final reordering and saving
prefixes = ['flt', params.coloc_prefix, params.eqtl_prefix, params.loop_prefix, 'gene']
flt_cols = [i for i in np.arange(master.shape[1]) if master.columns[i].startswith('flt.')]
gene_cols = [i for i in np.arange(master.shape[1]) if master.columns[i].startswith('gene.')]
coloc_cols = [i for i in np.arange(master.shape[1]) if master.columns[i].startswith('coloc.')]
eqtl_cols = [i for i in np.arange(master.shape[1]) if master.columns[i].startswith('eqtl.')]
main_cols = [i for i in np.arange(master.shape[1]) if master.columns[i].split('.')[0] not in prefixes]
loop_cols = [i for i in np.arange(master.shape[1]) if master.columns[i].startswith('loop.')]
final_master = master.iloc[:, main_cols + flt_cols + coloc_cols + loop_cols + eqtl_cols + gene_cols]

# remove a few 
#master.drop(['loop.id', 'coloc.id'], axis=1, errors='ignore', inplace=True)

# write out the master data
fn = os.path.join(params.outdir, 'master.tsv')
os.makedirs(params.outdir, exist_ok=True)
final_master.to_csv(fn, sep='\t', header=True, index=False)

#fn = os.path.join(params.outdir, 'master.xlsx')
#excel_master = master.sort_values('rsid').set_index('rsid')
#excel_master.to_excel(fn, na_rep='nan')
