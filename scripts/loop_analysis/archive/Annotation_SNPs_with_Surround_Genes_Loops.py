#!/usr/bin/env python
# coding: utf-8

# # Make the master table

import os 
import sys
import pybedtools as pbt
import pandas as pd
import numpy as np
import subprocess as sp
import json
os.chdir('/mnt/BioHome/jreyna/jreyna/projects/dchallenge/')
pbt.set_bedtools_path('/mnt/BioApps/bedtools/bin/')
bgzip = '/mnt/BioApps/tabix/tabix-0.2.6/bgzip'
tabix = '/mnt/BioApps/tabix/tabix-0.2.6/tabix'

bedpe_6cols = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB']
bedpe_10cols = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB', 'name', 'score', 'strand1', 'strand2']

## default values for the command line
#sys.argv = [0] * 8
#sys.argv[1] =  'results/main/2021_Nikhil_eQTL/Results/Colocalization/T1D_34012112_Gaulton/'
#sys.argv[1] += 'BLUEPRINT_eQTL_Monocyte/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed'
#sys.argv[2] = 'results/refs/ensembl/gencode.v19.annotation.bed'
#sys.argv[3] = 'results/main/2021_Nikhil_eQTL/Data/FitHiChIP_Loops/CM/FitHiChIP_L/FitHiChIP.interactions_FitHiC_Q0.01.bed'
#sys.argv[4] = 'results/refs/spp/SPP_D-Challenge_networks.xlsx'
#sys.argv[5] = 'results/refs/hg19/hg19.chrom.sizes'
#sys.argv[6] = 'results/main/2021_Nikhil_eQTL/Data/eqtl_sqtl_summ_stats/BLUEPRINT_eQTL/Monocyte.txt.gz'
#sys.argv[7] = 'results/main/loop_analysis/washU/'

# parsing the commandline arguments
coloc_fn = sys.argv[1]
genes_fn = sys.argv[2]
loop_fn = sys.argv[3]
spp_fn = sys.argv[4]
gs_fn = sys.argv[5]
eqtl_fn = sys.argv[6]
outdir = sys.argv[7]

# setting the output file names
os.makedirs(outdir, exist_ok=True)

# ## Load the colocalization data

# load the colocalization data
coloc = pd.read_table(coloc_fn)

# extract the most significant according the H4 
coloc_sig_df = coloc[coloc['pp_H4_Coloc_Summary'] > 0.75]
coloc_sig_df.rename(columns={'pos': 'end'}, inplace=True)
coloc_sig_df.loc[:, 'start'] = coloc_sig_df.loc[:, 'end'] - 1
coloc_sig_full = coloc_sig_df.copy(deep=True)

coloc_sig_df = coloc_sig_df[['chr', 'start', 'end', 'rs_id', 'variant_id']]
coloc_sig_df = coloc_sig_df.loc[~coloc_sig_df.duplicated(subset='rs_id'),]
coloc_sig_pbt = pbt.BedTool.from_dataframe(coloc_sig_df.iloc[:, 0:4]).sort()

#csnp_slop_pbt = coloc_sig_pbt.slop(b=500000, g=gs_fn)

# ## Load the gene data

# load the gencode coords
cols = ['chrom', 'start', 'end', 'strand', 'type', 'gene_id', 'gene_name']
gencode = pd.read_table(genes_fn, header=None, names=cols)

# extract just the genes
genes_df = gencode.loc[gencode.type.isin(['gene'])]
genes_df = genes_df.loc[~genes_df.duplicated(subset='gene_id'), :]
genes_df.loc[:, 'chrom'] = genes_df['chrom'].astype(str)
genes_df = genes_df.iloc[:, [0,1,2,6,5]]
genes_pbt = pbt.BedTool.from_dataframe(genes_df).sort()

# ## Find the closest gene

closest_gene = coloc_sig_pbt.closest(genes_pbt, d=True)
closest_gene = closest_gene.to_dataframe()
closest_gene = closest_gene.iloc[:, [3,7,8,9]]
closest_gene.columns = ['rs_id', 'cls_gname', 'cls_id', 'cls_dist']
closest_gene.head()

uniq_cls_gname = closest_gene.groupby(['rs_id']).cls_gname.apply(lambda x: ','.join(x))
uniq_cls_ids = closest_gene.groupby(['rs_id']).cls_id.apply(lambda x: ','.join(x))
uniq_cls_dist = closest_gene.groupby(['rs_id']).cls_dist.apply(lambda x: ','.join([str(i) for i in x]))
uniq_cls = pd.merge(uniq_cls_gname, uniq_cls_ids, left_index=True, right_index=True)
uniq_cls = pd.merge(uniq_cls, uniq_cls_dist, left_index=True, right_index=True)
uniq_cls.reset_index(inplace=True)


# In[9]:


uniq_cls


# ## Find all genes +/- 500kb

# In[10]:


# get a list of gene names within +- 500kb of the SNPs
fivekb_gnames = coloc_sig_pbt.slop(b=500000, g=gs_fn)
fivekb_gnames = fivekb_gnames.map(genes_pbt, c=4, o='collapse')
fivekb_gnames = fivekb_gnames.to_dataframe()
fivekb_gnames = fivekb_gnames.iloc[:, [3,4]]
fivekb_gnames.columns = ['rs_id', 'gene_name']

# get a list of gene ids within +- 500kb of the SNPs
fivekb_gids = coloc_sig_pbt.slop(b=500000, g=gs_fn)
fivekb_gids = fivekb_gids.map(genes_pbt, c=5, o='collapse')
fivekb_gids = fivekb_gids.to_dataframe()
fivekb_gids = fivekb_gids.iloc[:, [3,4]]
fivekb_gids.columns = ['rs_id', 'gene_id']

# merge the two above results
fivekb_genes = fivekb_gnames.merge(fivekb_gids, on='rs_id')
fivekb_genes.columns = ['rs_id', '5kb_gname', '5kb_gid']


# In[11]:


# get eQTL's
eqtls = pd.read_table(eqtl_fn)
eqtls.columns = ['eqtl_gname', 'nvar', 'shape1', 'shape2', 'dummy',
                 'sid', 'dist', 'npval', 'slope', 'ppval', 'bpval', 'qval']


# ## Get the loops

# In[12]:


# load the loop data
loops = pd.read_table(loop_fn)
tmp_loops = loops[['chr1', 's1', 'e1', 'chr2', 's2', 'e2']]
tmp_loops.rename(columns={'p': 'score'}, inplace=True)
tmp_loops.loc[:, 'name'] = '.'
tmp_loops.loc[:, 'score'] = loops['p']
tmp_loops.loc[:, 'strand1'] = '.'
tmp_loops.loc[:, 'strand2'] = '.'
loops = pbt.BedTool.from_dataframe(tmp_loops)
print('FitHiChIP found {} significant loops.'.format(tmp_loops.shape[0]))


# #### Find out SNP - 5kb gene pairs with loops

# In[13]:


fivekb_genesI = coloc_sig_pbt.slop(b=500000, g=gs_fn)
fivekb_genesI = fivekb_genesI.intersect(genes_pbt, wa=True, wb=True)
fivekb_genesI = fivekb_genesI.to_dataframe()
fivekb_genesI['start'] += 500000
fivekb_genesI['end'] -= 500000

# re-arranging to fit bedpe format
fivekb_genesI = fivekb_genesI.iloc[:, [0,1,2,4,5,6,3,7,8]]
fivekb_genesI['dummy'] = 'drop'

# loading into pbt
fivekb_genesI = pbt.BedTool.from_dataframe(fivekb_genesI)
fivekb_genesI = fivekb_genesI.pair_to_pair(loops, type='both',  **{'is':True})
fivekb_genesI = fivekb_genesI.to_dataframe(disable_auto_names=True, header=None)
fivekb_genesI = fivekb_genesI.iloc[:, [10, 11, 12, 13, 14, 15,6,7,8,17]]
fivekb_genesI.columns = bedpe_6cols + ['rs_id', 'gname', 'gid', 'fithic_qval']


# #### Find out SNP - eQTL gene pairs with loops

# ## Generate the master table

master = coloc_sig_full.copy()

# add sid which is the chr:position of the SNP
master['sid'] = master['chr'].str.replace('chr', '') + ':' +  master['end'].astype(str)

# add the closest gene
master = master.merge(uniq_cls, on='rs_id', how='left')

# add the +/- fivekb genes 
master = master.merge(fivekb_genes, on='rs_id', how='left')

# add the eQTL data
eqtl_genes = master.merge(eqtls[['sid', 'eqtl_gname']], on='sid')
eqtl_genes = eqtl_genes.groupby('rs_id').eqtl_gname.unique()
eqtl_genes = eqtl_genes.apply(lambda x: ','.join(x))
master = master.merge(eqtl_genes, on='rs_id', how='left')

new_colnames = [
 'rs_id',
 'variant_id',
 'chr',
 'start',
 'end',           
 'geneName',  
 'eqtl_gname',
 'cls_gname',
 '5kb_gname',
 'pp_H0_Coloc_Summary',
 'pp_H1_Coloc_Summary',
 'pp_H2_Coloc_Summary',
 'pp_H3_Coloc_Summary',
 'pp_H4_Coloc_Summary',           
 'dist',
 'pvalue',
 'FDR',
 'slope_snp',
 'ref',
 'alt',
 'AC',
 'AF',
 'AN',
 'slope_se_snp',
 'slope_gwas',
 'slope_se_gwas',
 'pval_nominal',
 'SampleSize',
 'sid',
 'cls_id',
 'cls_dist',
 '5kb_gid'
]
master = master.loc[:, new_colnames]
master.rename(columns={'geneName': 'coloc_gname',
                       'end': 'pos',
                       'eqtl_gname': 'eqtl_gnames', 
                       'cls_gname': 'cls_gnames', 
                       'cls_id': 'cls_ids'}, inplace=True)
master.drop(['start'], axis=1, inplace=True)

sg_with_loops = fivekb_genesI[['rs_id', 'gname']].values.tolist()
sg_with_loops = set([tuple(x) for x in sg_with_loops])

coloc_loops = []
eqtl_loops = []
closest_loops = []
fivekb_loops = []
for i, sr in master.iterrows():
    
    # analyze colocs
    gene = sr.coloc_gname
    bools = 0 
    if (sr.rs_id, gene) in sg_with_loops:
        bools = 1
    coloc_loops.append(bools)
    
    # analyze qtls
    bools = []
    for gene in sr.eqtl_gnames.split(','):
        eqtl_bool = 0 
        if (sr.rs_id, gene) in sg_with_loops:
            eqtl_bool = 1
        bools.append(eqtl_bool)
    eqtl_loops.append(bools)
    
    # analyze closest genes       
    bools = []
    for gene in sr.cls_gnames.split(','):
        cls_bool = 0 
        if (sr.rs_id, gene) in sg_with_loops:
            cls_bool = 1
        bools.append(cls_bool)
    closest_loops.append(bools)
    
    # analyze closest genes       
    bools = []
    for gene in sr['5kb_gname'].split(','):
        cls_bool = 0 
        if (sr.rs_id, gene) in sg_with_loops:
            cls_bool = 1
        bools.append(cls_bool)
    fivekb_loops.append(bools)

master['coloc_loops'] = coloc_loops
master['eqtl_loops'] = eqtl_loops
master['closest_loops'] = closest_loops
master['fivekb_loops'] = fivekb_loops
master['fivekb_loops'].iloc[2]

fn = os.path.join(outdir, 'master.tsv')
master.to_csv(fn, sep='\t', header=True, index=False)

fn = os.path.join(outdir, 'master.xlsx')
excel_master = master.sort_values('rs_id').set_index('rs_id')
excel_master.to_excel(fn)
