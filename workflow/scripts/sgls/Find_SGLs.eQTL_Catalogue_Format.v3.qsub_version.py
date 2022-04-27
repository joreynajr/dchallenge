#!/mnt/BioHome/jreyna/software/anaconda3/hichip-db/bin/ python

import os 
import sys
import pybedtools as pbt
import pandas as pd
import numpy as np
import subprocess as sp
import json
import argparse
pbt.helpers.set_tempdir('/mnt/BioHome/jreyna/tmp/')
pbt.set_bedtools_path('/mnt/BioApps/bedtools/bin/')
bgzip = '/mnt/BioApps/tabix/tabix-0.2.6/bgzip'
tabix = '/mnt/BioApps/tabix/tabix-0.2.6/tabix'

bedpe_6cols = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB']
bedpe_10cols = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB', 'name', 'score', 'strand1', 'strand2']

# Making a command line interface
parser = argparse.ArgumentParser()
parser.add_argument('--eqtl', type=str, required=True)
parser.add_argument('--coloc', type=str, required=True)
parser.add_argument('--loop', type=str, required=True)
parser.add_argument('--celltype', type=str, required=True)
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
params = parser.parse_args()

#print(params)

# # Parse the arguments
eqtl_fn = params.eqtl
coloc_fn = params.coloc
loop_fn = params.loop
celltype = params.celltype
gs_fn = params.gs
genes_fn = params.gene_ref
loop_slop = params.loop_slop
outdir = params.outdir

# # Load the colocalization data

coloc = pd.read_table(coloc_fn)
rename_coloc_cols = {
    'chr': 'coloc_chr',
    'pos': 'coloc_pos',
    'pp_H0_Coloc_Summary': 'ppH0',
    'pp_H1_Coloc_Summary': 'ppH1',
    'pp_H2_Coloc_Summary': 'ppH2',
    'pp_H3_Coloc_Summary': 'ppH3',
    'pp_H4_Coloc_Summary': 'ppH4',
    'rs_id': 'rsid',
    'variant_id': 'varid',
    'geneName': 'geneid',
    'dist': 'dist',
    'pvalue': 'eqtl_pval',
    'FDR': 'eqtl_fdr',
    'fdr': 'eqtl_fdr',
    'slope_snp': 'eqtl_slope',
    'ref': 'ref',
    'alt': 'alt',
    'AC': 'AC',
    'AF': 'AF',
    'AN': 'AN',
    'slope_se_snp': 'eqtl_slope_se',
    'slope_gwas': 'gwas_slope',
    'slope_se_gwas': 'gwas_slope_se',
    'pval_nominal': 'gwas_pval',
    'SampleSize': 'sample_size'}
coloc.rename(columns=rename_coloc_cols, inplace=True)

# extract the most significant according the H4 
coloc_sig_df = coloc.loc[~coloc.duplicated(subset=['rsid', 'geneid']),]

coloc_sig_full = coloc_sig_df.copy(deep=True)
coloc_sig_df.rename(columns={'pos': 'end'}, inplace=True)
coloc_sig_df.loc[:, 'start'] = coloc_sig_df.loc[:, 'end'] - 1

coloc_sig_df = coloc_sig_df[['chr', 'start', 'end', 'rsid', 'varid']]
coloc_sig_pbt = pbt.BedTool.from_dataframe(coloc_sig_df.iloc[:, 0:4]).sort()

print('There are {} colocalized SNP-gene pairs'.format(coloc_sig_df.shape[0]))

# # Load the gene data
print('# Load the gene data')

# load the gencode coords
cols = ['chrom', 'start', 'end', 'strand', 'type', 'geneid', 'genename']
gencode = pd.read_table(genes_fn, header=None, names=cols)

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

# # Find all genes +/- 500kb of a colocalized SNP
print('# Find all genes +/- 500kb of a colocalized SNP')

slop = 1000000
print("# Find all genes +/- {}kb".format(int(slop/1000)))

# get a list of gene names within +- 500kb of the SNPs
fivekb_genes = coloc_sig_pbt.slop(b=slop, g=gs_fn)
fivekb_genes = fivekb_genes.intersect(genes_pbt, wa=True, wb=True)
fivekb_genes = fivekb_genes.to_dataframe()
fivekb_genes = fivekb_genes.iloc[:, [0,1,2,4,5,6,3,7,8,9]]
fivekb_genes.columns = bedpe_6cols + ['rsid', 'genename', 'geneid', 'strand2']
fivekb_genes['strand1'] = '+'
fivekb_genes['name'] = fivekb_genes['rsid'] + '_' + fivekb_genes['genename']
fivekb_genes['score'] = '.'
new_order = bedpe_10cols + ['rsid', 'genename', 'geneid']
fivekb_genes = fivekb_genes[new_order]

fivekb_genes['startA'] += slop
fivekb_genes['endA'] -= slop
fivekb_genes['sid'] = fivekb_genes['chrA'].str.replace('chr', '') + ':' + fivekb_genes['endA'].astype(str)

print('There are {} colocalized snp-gene pairs within +/- 1000kb.'.format(fivekb_genes.shape[0]))


# # Find the closest gene
 
print('# Find the closest gene')

print('# Find the closest gene')
closest_gene = coloc_sig_pbt.closest(genes_pbt, d=True)
closest_gene = closest_gene.to_dataframe().iloc[:, [0,1,2,4,5,6,3,7,8,9]]
closest_gene.columns = bedpe_6cols + ['rsid', 'genename', 'geneid', 'dist']
closest_gene['sid'] = closest_gene['chrA'].str.replace('chr', '') + ':' + closest_gene['endA'].astype(str)
closest_gene.set_index(['sid', 'geneid'], inplace=True)

# # Load the H3K27ac HiChIP loops

print('# Load the H3K27ac HiChIP loops')

print('# Get the loops')

# load the loop data
loops = pd.read_table(loop_fn)
tmp_loops = loops[['chr1', 's1', 'e1', 'chr2', 's2', 'e2']]
tmp_loops.rename(columns={'p': 'score'}, inplace=True)
tmp_loops.loc[:, 'name'] = '.'
tmp_loops.loc[:, 'score'] = loops['p']
tmp_loops.loc[:, 'strand1'] = '.'
tmp_loops.loc[:, 'strand2'] = '.'
loops = pbt.BedTool.from_dataframe(tmp_loops)
print(tmp_loops.head())

# # Find SNP-Genes overlapping loops (both anchors + single anchor)

print('# Find SNP-Genes overlapping loops (both anchors + single anchor)')

# ### Find SNP-Gene overlapping a loop (both anchors)

print('### Find SNP-Gene overlapping a loop (both anchors)')

# re-arranging to fit bedpe format
fivekb_gloops = fivekb_genes.copy()

# loading into pbt
fivekb_gloops = pbt.BedTool.from_dataframe(fivekb_gloops)
fivekb_gloops = fivekb_gloops.pair_to_pair(loops, type='both', slop=loop_slop, **{'is':True})
fivekb_gloops = fivekb_gloops.to_dataframe(disable_auto_names=True, header=None)

if len(fivekb_gloops) > 0:
    fivekb_gloops_set = fivekb_gloops.iloc[:, [13,12]]
    fivekb_gloops_uniq = set([tuple(x) for x in fivekb_gloops_set.values.tolist()])
else:
    print('WARNING: Found no overlap between the genes and loops.')
    fivekb_gloops_uniq = set()

print('There are {} SNP-Gene pairs with a loop.'.format(len(fivekb_gloops_uniq)))


# ### Find coloc-SNP overlapping an anchor

print('### Find coloc-SNP overlapping an anchor')
colocSNP_anchors = loops.pairtobed(coloc_sig_pbt.slop(b=loop_slop, g=gs_fn), type='either')
colocSNP_anchors = colocSNP_anchors.to_dataframe(disable_auto_names=True, header=None)

if len(colocSNP_anchors) > 0:
    colocSNP_anchors_set = colocSNP_anchors.iloc[:, [10, 12]]
    colocSNP_anchors_set = [x[0].replace('chr', '') + ':' + str(x[1]) for x in colocSNP_anchors_set.values.tolist()]
else:
    print('WARNING: Found no overlap between coloc SNPs and loop anchors.')
    colocSNP_anchors_set = set()

print('There are {} SNPs which overlap a loop anchor.'.format(len(colocSNP_anchors_set)))

# # Construct master table

print('# Construct master table')

# begin making the master
master = fivekb_genes.copy()
master['sid'] = master['chrA'].str.replace('chr', '') + ':' + master['endA'].astype(str)

print('Master is starting with {} snp-gene pairs.'.format(master.shape[0]))

# #### Add eqtl results 
print('#### Add eqtl results')

# get SNP eQTL's; filtering required for eQTL Catalogue data
eqtls = pd.read_table(eqtl_fn, header=0)
if 'type' in eqtls.columns:
    eqtls = eqtls.loc[eqtls.type == 'SNP']

# renaming the main eqtl columns
eqtls_rename_cols = {
        eqtls.columns[params.eqtl_chr - 1]: 'chr',
        eqtls.columns[params.eqtl_pos - 1]: 'pos',
        eqtls.columns[params.eqtl_geneid - 1]: 'geneid',
        eqtls.columns[params.eqtl_beta - 1]: 'eqtl_beta',
        eqtls.columns[params.eqtl_pvalue - 1]: 'eqtl_pval',
        eqtls.columns[params.eqtl_dist - 1]: 'dist', 
        eqtls.columns[params.eqtl_fdr - 1]: 'eqtl_fdr'}
eqtls.rename(columns=eqtls_rename_cols, inplace=True)

# making an index list of main eqtl columns
main_eqtl_cols = [params.eqtl_chr, params.eqtl_pos, params.eqtl_geneid, 
                    params.eqtl_beta, params.eqtl_pvalue, params.eqtl_dist, params.eqtl_fdr]
main_eqtl_cols = [x - 1 for x in main_eqtl_cols] # converting to 0-based indexing

# making an index lists of extra eqtl columns
extra_eqtl_cols = list(range(eqtls.shape[1]))
for x in extra_eqtl_cols.copy():
    if x in main_eqtl_cols:
        extra_eqtl_cols.remove(x)

# reorder the eqtl dataset
eqtls = eqtls.iloc[:, main_eqtl_cols + extra_eqtl_cols]

# replace with new method of renaming above
#eqtls_rename_cols = {'chromosome': 'chr',
#    'position': 'pos',
#    'ref': 'ref',
#    'alt': 'alt',
#    'variant': 'variant',
#    'ma_samples': 'ma_samples',
#    'maf': 'maf',
#    'pvalue': 'eqtl_pval',
#    'pval': 'eqtl_pval',
#    'beta': 'eqtl_beta',
#    'se': 'eqtl_slop_se',
#    'ac': 'ac',
#    'an': 'an',
#    'r2': 'r2',
#    'gene_id': 'geneid',
#    'median_tpm': 'median_tpm',
#    'rsid': 'rsid',
#    'dist': 'dist',
#    'fdr': 'eqtl_fdr'}

eqtls['sid'] = eqtls['chr'].astype(str) + ':' + eqtls['pos'].astype(str)
eqtls['is_eqtl_pair'] = -1 # add column to filter on eqtl snp status
eqtls.loc[(eqtls['eqtl_fdr'] < 0.05), 'is_eqtl_pair'] = 1
eqtls.loc[(eqtls['eqtl_fdr'] > 0.05), 'is_eqtl_pair'] = 0

print('There are {} eQTLs.'.format(eqtls.shape[0]))

# need to use outer or else you exclude some eQTL's
master = master.merge(eqtls, on=['sid', 'geneid'], how='left')

# update eqtl snp status
master.loc[master.is_eqtl_pair.isna(), 'is_eqtl_pair' ] = 0 
master.loc[:, 'is_eqtl_pair'] = master.loc[:, 'is_eqtl_pair'].astype(int)

# add gene names to entries with a missing name (after adding eQTL info)
master.loc[master.genename.isna(), 'genename'] = master.loc[master.genename.isna(), 'geneid']

# add missing chrA, chrB, startA and endA data for the eQTL rows
master.loc[master.chrA.isna(), 'chrA'] = 'chr' + master.loc[master.chrA.isna(), 'sid'].str.replace(':[0-9]+', '')
master.loc[master.chrB.isna(), 'chrB'] = 'chr' + master.loc[master.chrB.isna(), 'sid'].str.replace(':[0-9]+', '')
master.loc[master.startA.isna(), 'startA'] = (master.loc[master.startA.isna(), 'sid'].str.replace('[0-9]+:', '')).astype(int)
master.loc[master.startA.isna(), 'startA'] -= 1     
master.loc[master.endA.isna(), 'endA'] = master.loc[master.endA.isna(), 'sid'].str.replace('[0-9]+:', '')

print('After outer merging with eqtls, master has {} snp-gene pairs.'.format(master.shape[0]))

# ### Add gene meta data 
print('### Add gene meta data')
#master.loc[~master.geneid.isna(), 'geneid'] = master.loc[~master.geneid.isna(), 'geneid']

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

master.rename(columns={'start': 'gene_start', 'end': 'gene_end', 'strand': 'gene_strand'}, inplace=True)

# ### Add info about closests gene

print('### Add info about closests gene')

# check for the closets gene
closets_check = [0] * master.shape[0]
for i, sr in master.iterrows():
    rs_gene = (sr.sid, sr.geneid)
    if rs_gene in closest_gene.index:
        closets_check[i] = 1
master['is_closest_gene'] = closets_check
print('Found {} closest genes.'.format(sum(master['is_closest_gene'])))

# ### Add colocalization data
print('### Add colocalization data')

# add colocalization data for SNP and is_coloc_snp columns
tmp_coloc = coloc_sig_full[['chr', 'pos',
 'ppH0',
 'ppH1',
 'ppH2',
 'ppH3',
 'ppH4',
 'rsid',
 'geneid',                       
 'ref',
 'alt',
 'AC',
 'AF',
 'AN',
 'gwas_slope',
 'gwas_slope_se',
 'gwas_pval']]
 ##'sample_size']]
# tmp_coloc.rename(columns={'slope_gwas': 'gwas_slope',
#                           'slope_se_gwas': 'gwas_slope_se',
#                           'pval_nominal': 'gwas_pval_nominal',
#                           'geneName': 'genename'}, inplace=True)
tmp_coloc['sid'] = tmp_coloc['chr'].str.replace('chr', '') + ':' +  tmp_coloc['pos'].astype(str) # NEW
master = master.merge(tmp_coloc, on=['sid', 'geneid'], how='left')


# ### Add column to filter on coloc snp status

print('### Add column to filter on coloc snp status')
master['is_coloc_pair'] = (~master['ppH4'].isna()).astype(int)

print('After left merging master with the colocalization table there are {} entries.'.format(master.shape[0]))
after_cnt = sum(master['is_coloc_pair'])
print('Checking if I have the correct number of colocalized SNPs: before {}; after {}.'.      format(coloc_sig_df.shape[0], after_cnt))

# ### Add loop data
print('### Add loop data')

# check for the loop gene
loop_check = [0] * master.shape[0]
for i, sr in master.iterrows():

    # check closest gene
    rs_gene = (sr.sid, sr.geneid)
    if rs_gene in fivekb_gloops_uniq:
        loop_check[i] = 1       

master['has_fithichip_loop'] = loop_check
print('There are {} SNP-Gene loops.'.format(sum(loop_check)))


# ### Add coloc-snp boolean column
print('Add coloc-snp boolean column')

master['has_colocSNP_anchor'] = master.sid.isin(colocSNP_anchors_set).astype(int)
print('There are {} coloc-SNP anchors.'.format(sum(master['has_colocSNP_anchor'])))

# #### Do the final reordering and saving

print("master.head()")
print(master.head())

# Include PC HiC data when it's available
#final_cols = ['sid', 
#    'rsid',
#    'geneid',
#    'genename',
#    'chrA',
#    'endA',
#    'startB',
#    'endB',
#    'is_eqtl_pair',
#    'is_coloc_pair',
#    'is_closest_gene',
#    'has_fithichip_loop',
#    'has_colocSNP_anchor',
#    'eqtl_pval', 
#    'eqtl_beta', 
#    'eqtl_fdr',
#    'dist',
#    'ppH0',
#    'ppH1',
#    'ppH2',
#    'ppH3',
#    'ppH4',
#    'ref_x', 
#    'alt_x',
#    'AC',
#    'AF',
#    'AN',
#    'gwas_slope',
#    'gwas_slope_se',
#    'gwas_pval',
#    'gene_start', 'gene_end', 'gene_strand']

final_cols = ['sid', 
    'geneid',
    'genename',
    'chrA',
    'endA',
    'startB',
    'endB',
    'is_eqtl_pair',
    'is_coloc_pair',
    'is_closest_gene',
    'has_fithichip_loop',
    'has_colocSNP_anchor',
    'eqtl_pval', 
    'eqtl_beta', 
    'eqtl_fdr',
    'dist',
    'ppH0',
    'ppH1',
    'ppH2',
    'ppH3',
    'ppH4',
    'AC',
    'AF',
    'AN',
    'gwas_slope',
    'gwas_slope_se',
    'gwas_pval',
    'gene_start', 'gene_end', 'gene_strand']


master = master[final_cols]
master.rename(columns={'chrA':'chrom',
                       'endA': 'snp_pos',
                       'ref_x': 'ref',
                       'alt_x': 'alt',
                       'startB': 'tss_start',
                       'endB': 'tss_end',
                       'genename': 'gene_name',
                       'geneid': 'geneid',
                       'pp_H0_Coloc_Summary': 'ppH0',
                       'pp_H1_Coloc_Summary': 'ppH1',
                       'pp_H2_Coloc_Summary': 'ppH2',
                       'pp_H3_Coloc_Summary': 'ppH3',
                       'pp_H4_Coloc_Summary': 'ppH4',
                       'pvalue': 'eqtl_pval',
                       'beta': 'eqtl_beta',
                       'fdr': 'eqtl_fdr',
                       'gwas_pval_nominal': 'gwas_pval'}, inplace=True)

master.sort_values(['chrom', 'snp_pos', 'tss_start', 'sid'], inplace=True)
master.snp_pos = master.snp_pos.astype(int)

# write out the master data
fn = os.path.join(outdir, 'master.tsv')
os.makedirs(outdir, exist_ok=True)
master.to_csv(fn, sep='\t', header=True, index=False)

#fn = os.path.join(outdir, 'master.xlsx')
#excel_master = master.sort_values('rsid').set_index('rsid')
#excel_master.to_excel(fn, na_rep='nan')
