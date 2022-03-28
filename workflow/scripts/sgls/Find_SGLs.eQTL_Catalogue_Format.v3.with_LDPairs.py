#!/usr/bin/env python
# coding: utf-8

import os 
import sys
import pybedtools as pbt
import pandas as pd
import numpy as np
import subprocess as sp
import json
from biolib import liftover as liftover_utility
from biolib import coordinates as biocoords
pbt.helpers.set_tempdir('/mnt/BioHome/jreyna/tmp/')
pbt.set_bedtools_path('/mnt/BioApps/bedtools/bin/')
bgzip = '/mnt/BioApps/tabix/tabix-0.2.6/bgzip'
tabix = '/mnt/BioApps/tabix/tabix-0.2.6/tabix'
res = 5000
pd.options.display.max_columns = 100


bedpe_6cols = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB']
bedpe_10cols = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB', 'name', 'score', 'strand1', 'strand2']


# # Parse the Snakemake arguments

# In[70]:


# snakemake.input = ['results/main/GRCh37/sgls/T1D_34012112_Gaulton/Quach_2016/monocyte_R848/eqtls.coloc_filtered.tsv.gz',
#               'results/main/coloc/eQTL_Catalogue/T1D_34012112_Gaulton/GRCh37/Quach_2016/monocyte_R848/',
#               'results/main/GRCh37/coloc/eQTL_Catalogue/T1D_34012112_Gaulton/Quach_2016/monocyte_R848/ldpairs/coloc_ld_snps.txt',
#               'results/main/h3k27ac_hichip/monocyte_naive/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed',
#               'results/refs/hg19/hg19.chrom.sizes',
#               'results/refs/gencode/v30/gencode.v30.annotation.grch37.bed']


# In[4]:


#eqtl_fn = snakemake.input[0]
#coloc_fn = snakemake.input[2]
#loop_fn = snakemake.input[3]
#celltype = snakemake.wildcards[3]
#gs_fn = snakemake.input[4]
#genes_fn = snakemake.input[5]
#loop_slop = snakemake.params[0]
#outdir = snakemake.output[0]

eqtl_fn = sys.argv[1] 
coloc_fn = sys.argv[2]
loop_fn = sys.argv[3]
celltype = sys.argv[4]
gs_fn = sys.argv[5]
genes_fn = sys.argv[6]
loop_slop = sys.argv[7]
outdir = sys.argv[8]

print("eqtl_fn:", eqtl_fn)
print("coloc_fn:", coloc_fn)
print("loop_fn:", loop_fn)
print("gs_fn:", gs_fn)
print("genes_fn:", genes_fn)
print("loop_slop:", loop_slop)
print("outdir:", outdir)

# # Load the colocalization data
# 

# load the colocalization data
coloc = pd.read_table(coloc_fn)

# rename the columns
rename_coloc_cols = {
    'chr': 'coloc_chr',
    'pos': 'coloc_pos',
    'pp_H0_Coloc_Summary': 'ppH0',
    'pp_H1_Coloc_Summary': 'ppH1',
    'pp_H2_Coloc_Summary': 'ppH2',
    'pp_H3_Coloc_Summary': 'ppH3',
    'pp_H4_Coloc_Summary': 'ppH4',
    'rs_id': 'coloc_rsid',
    'variant_id': 'varid',
    'geneName': 'geneid',
    'dist': 'dist',
    'pvalue': 'eqtl_pval',
    'FDR': 'eqtl_fdr',
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
    'SampleSize': 'sample_size',
    'rsID': 'coloc_rsid_dup',
    'ld_chr': 'chr',
    'ld_pos': 'pos',
    'ld_rsID': 'rsid'}
coloc.rename(columns=rename_coloc_cols, inplace=True)

# reorder the columns 
major_cols = ['chr', 'pos', 'rsid',
             'ppH0', 'ppH1', 'ppH2', 'ppH3', 'ppH4',
             'varid', 'geneid', 'dist',
             'eqtl_pval', 'eqtl_fdr', 'eqtl_slope',
             'ref', 'alt', 'AC', 'AF', 'AN',
             'eqtl_slope_se', 'gwas_slope', 'gwas_slope_se', 'gwas_pval',
             'coloc_chr', 'coloc_pos', 'coloc_rsid', 'LD']
coloc = coloc[major_cols]

# construct the sid
coloc.loc[:, 'sid'] = coloc['chr'].str.replace('chr', '') + ':' + coloc['pos'].astype(str) 


# In[14]:


# extract the most significant according the H4 
coloc_sig_df = coloc.loc[~coloc.duplicated(subset=['chr', 'pos', 'geneid']),]
coloc_sig_df.head()


# In[15]:


coloc_sig_df.shape


# In[16]:


# collect the chr and position of the main colocalized snps 
coloc_snps = coloc.loc[(coloc['chr'] == coloc['coloc_chr']) & (coloc['pos'] == coloc['coloc_pos'])] 
#coloc_snps = coloc_snps[['chr', 'pos']].values.tolist()
#coloc_snps = set([tuple(x) for x in coloc_snps])
coloc_snps = coloc_snps['sid'].unique()


# In[17]:


# make a pybedtool object for intersections
coloc_sig_full = coloc_sig_df.copy(deep=True)


# In[18]:


# get the bin start and end coordinates for each snp
coloc_bins = coloc_sig_full.pos.apply(biocoords.get_snp_bin, args=(res,))
coloc_starts, coloc_ends = list(zip(*coloc_bins.values))
coloc_sig_df.loc[:, 'start'] = coloc_starts
coloc_sig_df.loc[:, 'end'] = coloc_ends


# In[19]:


coloc_sig_df = coloc_sig_df[['chr', 'start', 'end', 'sid', 'rsid', 'varid']]
coloc_sig_pbt = pbt.BedTool.from_dataframe(coloc_sig_df.iloc[:, 0:5]).sort()

print('There are {} colocalized/LD SNP-gene pairs'.format(coloc_sig_df.shape[0]))


# # Load the gene data
# 

# In[20]:


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


# In[21]:


# get the bin start and end 
genes_pbt_pre = genes_df.copy()
gene_bins = genes_pbt_pre.start.apply(biocoords.get_snp_bin, args=(res,))
gene_starts, gene_ends = list(zip(*gene_bins.values))
genes_pbt_pre.loc[:, 'start'] = gene_starts 
genes_pbt_pre.loc[:, 'end'] = gene_ends


# In[22]:


# make a genes pbt for intersection
print("# make a genes pbt for intersection")
print(genes_df.head())
genes_pbt = pbt.BedTool.from_dataframe(genes_pbt_pre).sort()

print('\nThere are {} genes in this GTF-derived file.'.format(genes_df.shape[0]))


# In[23]:


genes_pbt_pre


# # Find all genes +/- 500kb of a colocalized/LD SNP

# In[24]:


print('# Find all genes +/- 500kb of a colocalized/LD SNP')


# In[25]:


slop = 1000000
print("# Find all genes +/- {}kb".format(int(slop/1000)))

# get a list of gene names within +- 500kb of the SNPs
fivekb_genes = coloc_sig_pbt.slop(b=slop, g=gs_fn)
fivekb_genes = fivekb_genes.intersect(genes_pbt, wa=True, wb=True)
fivekb_genes = fivekb_genes.to_dataframe()

fivekb_genes = fivekb_genes.iloc[:, [0,1,2,5,6,7,3,4,8,9,10]]
fivekb_genes.columns = bedpe_6cols + ['sid', 'rsid', 'genename', 'geneid', 'strand2']
fivekb_genes['strand1'] = '+'
fivekb_genes['name'] = fivekb_genes['rsid'] + '_' + fivekb_genes['genename']
fivekb_genes['score'] = '.'

# extract only major columns
major_cols = bedpe_10cols + ['sid', 'genename', 'geneid']
fivekb_genes = fivekb_genes[major_cols]


# In[26]:


print('There are {} colocalized snp-gene pairs within +/- 1000kb.'.format(fivekb_genes.shape[0]))


# # Find the closest gene
# 

# In[27]:


print('# Find the closest gene')
closest_gene = coloc_sig_pbt.closest(genes_pbt, d=True)
closest_gene = closest_gene.to_dataframe().iloc[:, [0,1,2,5,6,7,3,4,8,9,11]]
closest_gene.columns = bedpe_6cols + ['sid', 'rsid', 'genename', 'geneid', 'dist']
closest_gene.set_index(['sid', 'geneid'], inplace=True)


# # Load the H3K27ac HiChIP loops

# In[28]:


print('# Load the H3K27ac HiChIP loops')


# In[29]:


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


# # Find SG pairs overlapping loops (both anchors + single anchor)

# In[30]:


#print('# Find SNP-Genes overlapping loops (both anchors + single anchor)')


# ### Find SNP-Gene overlapping a loop (both anchors)

# In[31]:


print('### Find SNP-Gene overlapping a loop (both anchors)')


# In[32]:


# re-arranging to fit bedpe format
fivekb_gloops = fivekb_genes.copy()

# loading into pbt
fivekb_gloops = pbt.BedTool.from_dataframe(fivekb_gloops)
fivekb_gloops = fivekb_gloops.pair_to_pair(loops, type='both', slop=loop_slop, **{'is':True})
fivekb_gloops = fivekb_gloops.to_dataframe(disable_auto_names=True, header=None)


# In[33]:


fivekb_gloops


# In[34]:


if len(fivekb_gloops) > 0:
    fivekb_gloops_set = fivekb_gloops.iloc[:, [10, 12]]
    fivekb_gloops_uniq = set([tuple(x) for x in fivekb_gloops_set.values.tolist()])
else:
    print('WARNING: Found no overlap between the genes and loops.')
    fivekb_gloops_uniq = set()

print('There are {} SNP-Gene pairs with a loop.'.format(len(fivekb_gloops_uniq)))


# ### Find coloc-SNP overlapping an anchor

# In[35]:


print('### Find coloc-SNP overlapping an anchor')


# In[36]:


colocSNP_anchors = loops.pairtobed(coloc_sig_pbt.slop(b=loop_slop, g=gs_fn), type='either')
colocSNP_anchors = colocSNP_anchors.to_dataframe(disable_auto_names=True, header=None)


# In[37]:


if len(colocSNP_anchors) > 0:
    colocSNP_anchors_set = colocSNP_anchors.iloc[:, 13]
else:
    print('WARNING: Found no overlap between coloc SNPs and loop anchors.')
    colocSNP_anchors_set = set()

print('There are {} SNPs which overlap a loop anchor.'.format(len(colocSNP_anchors_set)))


# # Construct the master table

# In[59]:


print('# Construct master table')

# begin making the master
master = fivekb_genes.copy()

print('Master is starting with {} snp-gene pairs.'.format(master.shape[0]))

# adding back the rsid
coord_to_rsid = coloc[['sid', 'rsid']].set_index(['sid']).squeeze().to_dict()
master.loc[:, 'rsid'] = master.sid.apply(lambda x: coord_to_rsid[x])


# #### Add eqtl results 

# In[63]:


print("# #### Add eqtl results")

# get eQTL's
eqtls = pd.read_table(eqtl_fn, header=0)#, nrows=1000)
eqtls = eqtls.loc[eqtls.type == 'SNP']
eqtls.drop(['molecular_trait_object_id', 'molecular_trait_id', 'type'], axis=1, inplace=True)
eqtls_rename_cols = {'chromosome': 'chr',
    'position': 'pos',
    'ref': 'ref',
    'alt': 'alt',
    'variant': 'variant',
    'ma_samples': 'ma_samples',
    'maf': 'maf',
    'pvalue': 'eqtl_pval',
    'beta': 'eqtl_beta',
    'se': 'eqtl_slop_se',
    'ac': 'ac',
    'an': 'an',
    'r2': 'r2',
    'gene_id': 'geneid',
    'median_tpm': 'median_tpm',
    'rsid': 'rsid',
    'dist': 'dist',
    'fdr': 'eqtl_fdr'}
eqtls.rename(columns=eqtls_rename_cols, inplace=True)
#eqtls['pos'] = eqtls['pos'].astype(int)
eqtls['sid'] = eqtls['chr'].astype(str) + ':' + eqtls['pos'].astype(str).str.replace('.0', '')
eqtls['is_eqtl_pair'] = -1 # add column to filter on eqtl snp status
eqtls.loc[(eqtls['eqtl_fdr'] <= 0.05), 'is_eqtl_pair'] = 1
eqtls.loc[(eqtls['eqtl_fdr'] > 0.05), 'is_eqtl_pair'] = 0


# In[65]:


print('There are {} eQTLs.'.format(eqtls.shape[0]))

# need to use outer or else you exclude some eQTL's
master = master.merge(eqtls, on=['sid', 'geneid'], how='left')


# In[ ]:


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


# ### Add the original gene coordinates 

# In[42]:


print('### Add the original gene coordinates ')


# In[43]:


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


# In[44]:


master


# ### Add info about closests gene

# In[45]:


print('### Add info about closests gene')


# In[46]:


print("# #### Add info about closests gene")

# check for the closets gene
closets_check = [0] * master.shape[0]
for i, sr in master.iterrows():
    rs_gene = (sr.sid, sr.geneid)
    if rs_gene in closest_gene.index:
        closets_check[i] = 1
master['is_closest_gene'] = closets_check
print('Found {} closest genes.'.format(sum(master['is_closest_gene'])))


# ### Add colocalization data

# In[47]:


print('### Add colocalization data')


# In[48]:


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
tmp_coloc['sid'] = tmp_coloc['chr'].str.replace('chr', '') + ':' +  tmp_coloc['pos'].astype(str) # NEW
master = master.merge(tmp_coloc, on=['sid', 'geneid'], how='left')


# In[49]:


master


# ### Add column to filter on coloc snp status

# In[50]:


print('### Add column to filter on coloc snp status')


# In[51]:


master['is_coloc_pair'] = (~master['ppH4'].isna()).astype(int)

print('After left merging master with the colocalization table there are {} entries.'.format(master.shape[0]))
after_cnt = sum(master['is_coloc_pair'])
print('Checking if I have the correct number of colocalized SNPs: before {}; after {}.'.      format(coloc_sig_df.shape[0], after_cnt))


# ### Add loop data

# In[52]:


print('# #### Add loop data')

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

# In[53]:


print('Add coloc-snp boolean column')


# In[54]:


master['has_colocSNP_anchor'] = master.sid.isin(colocSNP_anchors_set).astype(int)
print('There are {} coloc-SNP anchors.'.format(sum(master['has_colocSNP_anchor'])))


# ### Add LD SNP boolean column

# In[55]:


# add a True/False (1/0) column for the presence of a direct coloc snp (non-LD)
is_ld_snps = [x not in coloc_snps for x in master.sid]
master['is_ld_snp'] = [int(x) for x in is_ld_snps]


# ### Final ordering 

# In[57]:


# #### Do the final reordering and saving

# Include PC HiC data when it's available
final_cols = ['sid', 
    'rsid',
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
    'is_ld_snp', 
    'eqtl_pval', 
    'eqtl_beta', 
     #'eqtl_fdr',
    'dist',
    'ppH0',
    'ppH1',
    'ppH2',
    'ppH3',
    'ppH4',
    'gene_start',
    'gene_end',
    'ref_x', 
    'alt_x',
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

master.sort_values(['chrom', 'snp_pos', 'tss_start', 'rsid'], inplace=True)
master.snp_pos = master.snp_pos.astype(int)

# write out the master data
fn = os.path.join(outdir, 'master.tsv')
os.makedirs(outdir, exist_ok=True)
master.to_csv(fn, sep='\t', header=True, index=False)

#fn = os.path.join(outdir, 'master.xlsx')
#excel_master = master.sort_values('rsid').set_index('rsid')
#excel_master.to_excel(fn, na_rep='nan')


# In[58]:


master


# In[ ]:





# In[ ]:




#master[master.is_eqtl_pair == 1].head()master[master.is_coloc_pair == 1].head()master[master.has_colocSNP_anchor == 1].head()master[master.has_fithichip_loop == 1].head()#########################################################################################
## Make WashU files 
##########################################################################################
#
#def bedpe_to_WashU_longrange(fn, df):
#    """
#       Convert from a loop bedpe file into WashU longrange, 
#       includes bgzip and tabix of the fn. 
#       
#       Params
#       -------
#       fn: str
#           path to the longrange output file (without gz)
#           
#       df: dataframe
#           columns 1-6 are as expected and column 7 is the p or q-value. 
#           
#       Output
#       ------
#       gzfn: str
#           path to the longrange with bgzip compression
#       tabix_fn: str
#           path to the index of the longrange file
#           
#   """
#
#    # parsing the data into WashU longrage format
#    data = []
#    for sr in df.values.tolist():
#
#        # calculate the -log(FDR)
#        qval = -np.log(sr[6])
#
#        # get the first pair data
#        second_pair_str = '{}:{}-{},{:.5f}'.format(*sr[3:6], qval)
#        first_row = sr[0:3] + [second_pair_str]
#
#        # get the second pair data
#        first_pair_str = '{}:{}-{},{:.5f}'.format(*sr[0:3], qval)
#        second_row = sr[3:6] + [first_pair_str]
#
#        # add each data row
#        data.append(first_row)
#        data.append(second_row)
#
#    data = sorted(data, key=lambda x: (x[0], x[1], x[2]))
#
#    # writing out the data
#    with open(fn, 'w') as f:
#        for line in data:
#            info = [str(x) for x in line]
#            info = '\t'.join(info)
#            f.write(info + '\n')
#           
#    # run bgzip
#    cmd = '{} {}'.format(bgzip, fn)
#    print(cmd)
#    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)
#
#    out, err = job.communicate()
#    print('out:', out.decode())
#    print('err:', err.decode())
#
#    # run tabix
#    lrange_gzfn = fn + '.gz'
#    cmd = '{} -f {}'.format(tabix, lrange_gzfn)
#    print(cmd)
#    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)
#
#    out, err = job.communicate()
#    print('out:', out.decode())
#    print('err:', err.decode())
#
#    print('Created the gzfn: {}'.format(fn + '.gz'))
#    print('Created the tabix: {}'.format(fn + '.gz.tbi'))
#
#def bed_WashU_bedgz(fn, df):
#    """
#       Convert from a bed dataframe into WashU longrange file 
#       includes bgzip and tabix of the fn. 
#       
#       Params
#       -------
#       fn: str
#           path to the longrange output file (without gz)
#           
#       df: dataframe
#           columns 1-3 are as expected and column 7 is the p or q-value. 
#           
#       Output
#       ------
#       gzfn: str
#           path to the longrange with bgzip compression
#       tabix_fn: str
#           path to the index of the longrange file
#           
#   """
#
#    # parsing the data into WashU longrage format
#    data = []
#    for sr in df.values.tolist():
#        data.append(sr[0:4])
#    data = sorted(data, key=lambda x: (x[0], x[1], x[2]))
#
#    # writing out the data
#    with open(fn, 'w') as f:
#        for line in data:
#            info = [str(x) for x in line]
#            info = '\t'.join(info)
#            f.write(info + '\n')
#
#    # run bgzip
#    cmd = '{} {}'.format(bgzip, fn)
#    print(cmd)
#    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)
#
#    out, err = job.communicate()
#    print('out:', out.decode())
#    print('err:', err.decode())
#
#    # run tabix
#    gzfn = fn + '.gz'
#    cmd = '{} -f {}'.format(tabix, gzfn)
#    print(cmd)
#    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)
#
#    out, err = job.communicate()
#    print('out:', out.decode())
#    print('err:', err.decode())
#
#    print('Created the gzfn: {}'.format(fn + '.gz'))
#    print('Created the tabix: {}'.format(fn + '.gz.tbi'))
#
#def bed_to_WashU_refbed(fn, df):
#    """ 
#       Convert from a bed dataframe into WashU longrange file 
#       includes bgzip and tabix of the fn. 
#       
#       Params
#       -------
#       fn: str
#           path to the longrange output file (without gz)
#           
#       df: dataframe
#           columns 1-3 are as expected and column 7 is the p or q-value. 
#           
#       Output
#       ------
#       gzfn: str
#           path to the longrange with bgzip compression
#       tabix_fn: str
#           path to the index of the longrange file
#           
#   """
#
#    # parsing the data into WashU longrage format
#    data = df.values.tolist()
#    data = sorted(data, key=lambda x: (x[0], x[1], x[2]))
#
#    # writing out the data
#    with open(fn, 'w') as f:
#        for line in data:
#            info = [str(x) for x in line]
#            info = '\t'.join(info)
#            f.write(info + '\n')
#
#    # run bgzip
#    cmd = '{} -f {}'.format(bgzip, fn) 
#    print(cmd)
#    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)
#
#    out, err = job.communicate()
#    print('out:', out.decode())
#    print('err:', err.decode())
#
#    # run tabix
#    gzfn = fn + '.gz'
#    cmd = '{} {}'.format(tabix, gzfn)
#    print(cmd)
#    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)
#
#    out, err = job.communicate()
#    print('out:', out.decode())
#    print('err:', err.decode())
#
#    print('Created the gzfn: {}'.format(fn + '.gz'))
#    print('Created the tabix: {}'.format(fn + '.gz.tbi'))# make the refbed link for genes (status: running)
#final_sg_cols = ['chrom', 'gene_start' ,'gene_end', 'gene_name', 'gene_strand']
#final_sg_genes = master.loc[(master.has_fithichip_loop == 1), final_sg_cols]
#
#final_sg_genes.gene_start = final_sg_genes.gene_start.astype(int)
#final_sg_genes.gene_end = final_sg_genes.gene_end.astype(int)
#
#final_sg_genes = final_sg_genes.loc[~final_sg_genes.duplicated()]
#final_sg_genes['chr'] = final_sg_genes['chrom'] 
#final_sg_genes['transcript_start'] = final_sg_genes['gene_start']
#final_sg_genes['transcript_stop'] = final_sg_genes['gene_end']
#final_sg_genes['translation_start'] = final_sg_genes['gene_start']
#final_sg_genes['translation_stop'] = final_sg_genes['gene_end']
#final_sg_genes['strand'] = final_sg_genes['gene_strand']
#final_sg_genes['gene_name'] = final_sg_genes['gene_name']
#final_sg_genes['transcript_id'] = final_sg_genes['gene_name']
#final_sg_genes['type'] = 'coding'
#final_sg_genes['exon_gene_start'] = final_sg_genes['gene_start']
#final_sg_genes['exon_stops'] = final_sg_genes['gene_end']
#refcols = ['chr', 'transcript_start', 'transcript_stop', 'translation_start',
#          'translation_stop', 'strand', 'gene_name', 'transcript_id',
#          'type', 'exon_gene_start', 'exon_stops']
#final_sg_genes = final_sg_genes.loc[:, refcols]
#sg_genes_fn = os.path.join(outdir, 'gs_genes.bed')
#bed_to_WashU_refbed(sg_genes_fn, final_sg_genes)
#
#
## make the longrange link for raw fithichip data (status: running)
#loop_gz = os.path.abspath(loop_fn.replace('.bed', '_WashU.bed.gz'))
#loop_tbi = os.path.abspath(loop_fn.replace('.bed', '_WashU.bed.gz.tbi'))
#loop_gz_link = os.path.join(outdir, os.path.basename(loop_gz))
#loop_tbi_link = os.path.join(outdir, os.path.basename(loop_tbi))
#if not os.path.exists(loop_gz_link):
#    os.link(loop_gz, loop_gz_link)
#    os.link(loop_tbi, loop_tbi_link)# ## All 5kb Washu Files
#
## make the longrange link for snp-gene fivekb pairs (status: running)
##fivekb_lrange = fivekb_lrange[fivekb_lrange.rsid.notna()].reset_index(drop=True)
#fivekb_lrange = fivekb_genes.copy()
#
## convert full for viz
#fivekb_lrange = fivekb_lrange[['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB']]
#fivekb_lrange.iloc[:, 1] -= 1 
#fivekb_lrange['score'] = 0.01
#
#fivekb_snp_gene_pairs_fn = os.path.join(outdir, '5kb.snp_gene_pairs.bed')
#bedpe_to_WashU_longrange(fivekb_snp_gene_pairs_fn, fivekb_lrange)
#
## plot all snp-gene pairs with a loop
#fivekb_gloops_lrange = fivekb_gloops.iloc[:, [11,12,13,14,15,16,18]]
#
#fivekb_snp_gene_loops_fn = os.path.join(outdir, '5kb.snp_gene_loops.bed')
#bedpe_to_WashU_longrange(fivekb_snp_gene_loops_fn, fivekb_gloops_lrange)# ## eQTL WashU Files
#
## make the longrange link for snp-gene eQTL pairs (status: running)
#eqtl_lrange = master.loc[master.is_eqtl_pair == 1]
#eqtl_lrange = eqtl_lrange[eqtl_lrange.rsid.notna()].reset_index(drop=True)
#
## convert full for viz
#eqtl_snp_gene_pairs = eqtl_lrange[['chrom', 'snp_pos', 'snp_pos', 'chrom', 'tss_start', 'tss_end', 'genename']]
#eqtl_snp_gene_pairs.iloc[:, 1] -= 1 
#eqtl_snp_gene_pairs['score'] = 0.01
#
#eqtl_snp_gene_pairs_fn = os.path.join(outdir, 'eqtl.snp_gene_pairs.bed')
#bedpe_to_WashU_longrange(eqtl_snp_gene_pairs_fn, eqtl_snp_gene_pairs.iloc[:, [0,1,2,3,4,5,7]])
#
## make the bed for egenes only (status: running)
#eqtl_genes = eqtl_lrange.loc[:, ['chrom', 'gene_start', 'gene_end', 'gene_name']]
#eqtl_genes_only_fn = os.path.join(outdir, 'eqtl.genes_only.bed')
#bed_WashU_bedgz(eqtl_genes_only_fn, eqtl_genes)# make the bed for eSNPs only (status: running)
#eqtl_snps = eqtl_lrange.loc[:, ['chrom', 'snp_pos', 'snp_pos', 'rsid']]
#eqtl_snps.iloc[:, 1] -= 1
#eqtl_snps = pbt.BedTool.from_dataframe(eqtl_snps)
#eqtl_snps = eqtl_snps.slop(b=500, g=gs_fn).to_dataframe()
#eqtl_snps_only_fn = os.path.join(outdir, 'eqtl.snps_only.bed')
#bed_WashU_bedgz(eqtl_snps_only_fn, eqtl_snps)
#
## # make the longrange link for snp-gene eQTL loops (status: didn't add loop coordinates to master.)
## eqtl_loops = master.loc[master.is_eqtl_pair == 1]# ## Coloc WashU Files
#
## make the bed for colocalized SNPs (status: running)
#final_snps = coloc_sig_df.copy()
#final_snps = final_snps.loc[~final_snps.duplicated()]
#final_snps = pbt.BedTool.from_dataframe(final_snps)
#final_snps = final_snps.slop(b=500, g=gs_fn)
#final_snps = final_snps.to_dataframe()[0:4]
#
#coloc_snps_only_fn = os.path.join(outdir, 'coloc.snps_only.bed')
#bed_WashU_bedgz(coloc_snps_only_fn, final_snps)
#
## make the bed for colocalized genes (status: running)
#final_coloc_snp_genes_pairs = master[(master.is_coloc_pair == 1)]
#final_coloc_genes = final_coloc_snp_genes_pairs[['chrom', 'gene_start' ,'gene_end', 'gene_name']]
#coloc_genes_only_fn = os.path.join(outdir, 'coloc.genes_only.bed')
#bed_WashU_bedgz(coloc_genes_only_fn, final_coloc_genes)# make the longrange for colocalized snp-gene pairs (status: running)
#final_coloc_snp_genes_pairs_out = final_coloc_snp_genes_pairs[['chrom', 'snp_pos', 'snp_pos', 'chrom', 
#                                                              'tss_start', 'tss_end', 'gene_name']]
#final_coloc_snp_genes_pairs_out.iloc[:, 1] -= 1
#final_coloc_snp_genes_pairs_out.iloc[:, 6] = 0.01
#
#coloc_snp_gene_pairs_fn = os.path.join(outdir, 'coloc.snp_gene_pairs.bed')
#bedpe_to_WashU_longrange(coloc_snp_gene_pairs_fn, final_coloc_snp_genes_pairs_out)# make the longrange for colocalized snp-gene loops (status: running)
#final_coloc_snp_genes_loops = fivekb_gloops[(fivekb_gloops[10].isin(final_coloc_snp_genes_pairs['rsid'])) & 
#             (fivekb_gloops[12].isin(final_coloc_snp_genes_pairs['geneid']))]
#final_coloc_snp_genes_loops_out = final_coloc_snp_genes_loops.iloc[:, [14,15,16,17,18,19,21]]
#
#coloc_snp_gene_loops_fn = os.path.join(outdir, 'coloc.snp_gene_loops.bed')
#bedpe_to_WashU_longrange(coloc_snp_gene_loops_fn, final_coloc_snp_genes_loops_out)
#
#final_coloc_snp_genes_loops_out# make the loop anchors as bed files (status: developing)
#left = final_coloc_snp_genes_loops[[14,15,16,6]].T.reset_index(drop=True).T
#left[3] = 'L-' + left[3] 
#right = final_coloc_snp_genes_loops[[17,18,19,6]].T.reset_index(drop=True).T
#right[3] = 'R-' + right[3] 
#anchors = pd.concat([left, right], ignore_index=True, axis=0)coloc_anchors_fn = os.path.join(outdir, 'coloc.anchors_only.bed')
#bed_WashU_bedgz(coloc_anchors_fn, anchors)
#
#
## # make the hub json file
#
##gwas, cline = coloc_fn.split('/')[5:7]
#
#print("# make the hub json file")
#
#gene_refbed_json = {'type': 'refbed',
#                   'filename': os.path.basename(sg_genes_fn) + '.gz',
#                   'name': 'Gencode V19',
#                   'showOnHubLoad': True
#                  }
#
#orig_loops_json = {'type': 'longrange',
#                  'filename': os.path.basename(loop_gz_link),
#                  'name': 'Original Loops',
#                  'options': {'displayMode': 'arc', 'color':'red'},
#                  'showOnHubLoad': True
#               }
#
## fivekb_snp_gene_pairs_json = {'type': 'longrange',
##                    'filename': os.path.basename(fivekb_snp_gene_pairs_fn) + '.gz',
##                    'name': '5kb SNP-Gene Pairs',
##                    'options': {'displayMode': 'arc', 'color':'purple', 'height': 200},
##                    'showOnHubLoad': False
##                 }
#
## fivekb_snp_gene_loops_json = {'type': 'longrange',
##                    'filename': os.path.basename(fivekb_snp_gene_loops_fn) + '.gz',
##                    'name': '5kb SNP-Gene Loops',
##                    'options': {'displayMode': 'arc', 'color':'red', 'height': 200},
##                    'showOnHubLoad': False
##                 }
#
#
#eqtl_snp_gene_pairs_json = {'type': 'longrange',
#                'filename': os.path.basename(eqtl_snp_gene_pairs_fn) + '.gz',
#                'name': 'eQTL SNP-Gene Pairs',
#                'options': {'displayMode': 'arc', 'color':'purple', 'height': 200},
#                'showOnHubLoad': True
#               }
#
#
#eqtl_snps_only_json = {'type': 'bed',
#                'filename': os.path.basename(eqtl_snps_only_fn) + '.gz',
#                'name': 'eQTL SNPs only',
#                'options': {'color':'purple'},
#                'showOnHubLoad': True
#
#               }
#
#
#eqtl_genes_only_json = {'type': 'bed',
#                'filename': os.path.basename(eqtl_genes_only_fn) + '.gz',
#                'name': 'eQTL genes only',
#                'options': {'color':'purple'},
#                'showOnHubLoad': True
#               }
#
#
#coloc_snp_gene_pairs_json = {'type': 'longrange',
#                'filename': os.path.basename(coloc_snp_gene_pairs_fn) + '.gz',
#                'name': 'coloc SNP-Gene Pairs',
#                'options': {'displayMode': 'arc', 'color':'purple', 'height': 200},
#                'showOnHubLoad': True
#               }
#
#coloc_snp_gene_loops_json = {'type': 'longrange',
#                'filename': os.path.basename(coloc_snp_gene_loops_fn) + '.gz',
#                'name': 'coloc SNP-Gene Loops',
#                'options': {'displayMode': 'arc', 'color':'red', 'height': 200},
#                'showOnHubLoad': True
#               }
#
#coloc_snps_only_json = {'type': 'bed',
#                'filename': os.path.basename(coloc_snps_only_fn) + '.gz',
#                'name': 'coloc SNPs only',
#                'options': {'color':'purple'},
#                'showOnHubLoad': False
#               }
#
#
#coloc_genes_only_json = {'type': 'bed',
#                'filename': os.path.basename(coloc_genes_only_fn) + '.gz',
#                'name': 'coloc genes only',
#                'options': {'color':'purple'},
#                'showOnHubLoad': False
#               }
#
#
#coloc_anchors_only_json = {'type': 'bed',
#                'filename': os.path.basename(coloc_anchors_fn) + '.gz',
#                'name': 'coloc loop anchors only',
#                'options': {'color':'red'},
#                'showOnHubLoad': False
#               }
#
#
##hub_json = [orig_loops_json, sg_pairs_json, sg_loops_json, sg_snps_json, sg_genes_json]
#hub_json = [gene_refbed_json, 
#           #fivekb_snp_gene_pairs_json,
#           #fivekb_snp_gene_loops_json,
#           eqtl_snps_only_json,
#           eqtl_genes_only_json,
#           eqtl_snp_gene_pairs_json,
#           coloc_snps_only_json,
#           coloc_genes_only_json,
#           coloc_snp_gene_pairs_json, 
#           coloc_anchors_only_json,
#           coloc_snp_gene_loops_json,
#           orig_loops_json, 
#          ]
#
#hub_json_fn = os.path.join(outdir, 'hub.config.json')
#with open(hub_json_fn, 'w') as f:
#    f.write(json.dumps(hub_json, indent=4))
## In[ ]:
#
#
#
#
#
## In[ ]:
#
#
#
#
