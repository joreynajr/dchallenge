#!/usr/bin/env python
# coding: utf-8

# %load ./scripts/loop_analysis/Annotate_Coloc_SNPs_with_Nearby_Genes_Loops.py
#!/usr/bin/env python

# # Make the master table
import os 
import sys
import pybedtools as pbt
import pandas as pd
import numpy as np
import subprocess as sp
import json
pbt.helpers.set_tempdir('/mnt/BioHome/jreyna/tmp/')
pbt.set_bedtools_path('/mnt/BioApps/bedtools/bin/')
bgzip = '/mnt/BioApps/tabix/tabix-0.2.6/bgzip'
tabix = '/mnt/BioApps/tabix/tabix-0.2.6/tabix'

bedpe_6cols = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB']
bedpe_10cols = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB', 'name', 'score', 'strand1', 'strand2']

# ## default values for the command line
#sys.argv = [0] * 100
#sys.argv[1] =  'results/main/2021_Nikhil_eQTL/Results/Colocalization/T1D_34012112_Gaulton/'
#sys.argv[1] += 'DICE_eQTL_CD4_NAIVE/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed'
#sys.argv[2] = 'results/refs/ensembl/gencode.v19.annotation.bed'
#sys.argv[3] = 'results/main/2021_Nikhil_eQTL/Data/FitHiChIP_Loops/CD4N/FitHiChIP_L/FitHiChIP.interactions_FitHiC_Q0.01.bed'
#sys.argv[5] = 'results/refs/hg19/hg19.chrom.sizes'
#sys.argv[6] = 'results/main/2021_Nikhil_eQTL/Data/eqtl_sqtl_summ_stats/DICE_eQTL/CD4_NAIVE.txt.gz'
#sys.argv[7] = 'results/main/pc_hic/data/PCHiC_peak_matrix_cutoff5.tsv'
#sys.argv[8] = 'CD4_NAIVE'
#sys.argv[9] = 'results/main/loop_analysis/test/'

#########################################################################################
# parsing the commandline arguments
#########################################################################################

eqtl_fn = sys.argv[1]
coloc_fn = sys.argv[2]
loop_fn = sys.argv[3]
celltype = sys.argv[4]
gs_fn = sys.argv[5]
genes_fn = sys.argv[6]
outdir = sys.argv[7]

#########################################################################################
# Load the colocalization data
#########################################################################################

# load the colocalization data
coloc = pd.read_table(coloc_fn)

# extract the most significant according the H4 
coloc_sig_df = coloc[coloc['pp_H4_Coloc_Summary'] > 0.75]
coloc_sig_df = coloc_sig_df.loc[~coloc_sig_df.duplicated(subset=['rs_id', 'geneName']),]

coloc_sig_full = coloc_sig_df.copy(deep=True)
coloc_sig_df.rename(columns={'pos': 'end'}, inplace=True)
coloc_sig_df.loc[:, 'start'] = coloc_sig_df.loc[:, 'end'] - 1

coloc_sig_df = coloc_sig_df[['chr', 'start', 'end', 'rs_id', 'variant_id']]
coloc_sig_pbt = pbt.BedTool.from_dataframe(coloc_sig_df.iloc[:, 0:4]).sort()

print('There are {} colocalized SNP-gene pairs'.format(coloc_sig_df.shape[0]))

#########################################################################################
# ## Load the gene data
#########################################################################################
print('# Load the gene data')

# load the gencode coords
cols = ['chrom', 'start', 'end', 'strand', 'type', 'gene_id', 'gname']
gencode = pd.read_table(genes_fn, header=None, names=cols)

# extract just the genes
genes_df = gencode.loc[gencode.type.isin(['gene'])]
genes_df = genes_df.loc[~genes_df.duplicated(subset='gene_id'), :]
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

print('There are {} genes in this GTF-derived file.'.format(genes_df.shape[0]))

#########################################################################################
# Find all genes +/- 500kb
#########################################################################################
print("# Find all genes +/- 500kb")

# get a list of gene names within +- 500kb of the SNPs
fivekb_genes = coloc_sig_pbt.slop(b=500000, g=gs_fn)
fivekb_genes = fivekb_genes.intersect(genes_pbt, wa=True, wb=True)
fivekb_genes = fivekb_genes.to_dataframe()
fivekb_genes = fivekb_genes.iloc[:, [0,1,2,4,5,6,3,7,8,9]]

fivekb_genes.columns = bedpe_6cols + ['rs_id', 'gname', 'gid', 'strand2']
fivekb_genes['strand1'] = '+'
fivekb_genes['name'] = fivekb_genes['rs_id'] + '_' + fivekb_genes['gname']
fivekb_genes['score'] = '.'
new_order = bedpe_10cols + ['rs_id', 'gname', 'gid']
fivekb_genes = fivekb_genes[new_order]

fivekb_genes['startA'] += 500000
fivekb_genes['endA'] -= 500000
fivekb_genes['sid'] = fivekb_genes['chrA'].str.replace('chr', '') + ':' + fivekb_genes['endA'].astype(str)

print('There are {} colocalized snp-gene pairs within +/- 5kb.'.format(fivekb_genes.shape[0]))

#########################################################################################
# Find the closest gene
#########################################################################################

closest_gene = coloc_sig_pbt.closest(genes_pbt, d=True)
closest_gene = closest_gene.to_dataframe().iloc[:, [0,1,2,4,5,6,3,7,8,9]]
closest_gene.columns = bedpe_6cols + ['rs_id', 'gname', 'gid', 'dist']
closest_gene['sid'] = closest_gene['chrA'].str.replace('chr', '') + ':' + closest_gene['endA'].astype(str)
closest_gene.set_index(['sid', 'gname'], inplace=True)


#########################################################################################
# Get the loops
#########################################################################################

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


# ## Load Promoter Capture HiC Data

# re-arranging to fit bedpe format
fivekb_gloops = fivekb_genes.copy()

# loading into pbt
fivekb_gloops = pbt.BedTool.from_dataframe(fivekb_gloops)
fivekb_gloops = fivekb_gloops.pair_to_pair(loops, type='both', slop=7500, **{'is':True})
fivekb_gloops = fivekb_gloops.to_dataframe(disable_auto_names=True, header=None)

fivekb_gloops_set = fivekb_gloops.iloc[:, [13,11]]
fivekb_gloops_uniq = set([tuple(x) for x in fivekb_gloops_set.values.tolist()])

print('There are {} SNP-Gene pairs with a loop.'.format(len(fivekb_gloops_uniq)))


#########################################################################################
# Construct master table
#########################################################################################

# begin making the master
master = fivekb_genes.copy()
master['sid'] = master['chrA'].str.replace('chr', '') + ':' + master['endA'].astype(str)

print('Master is starting with {} snp-gene pairs.'.format(master.shape[0]))

# #### Add eqtl results 
print("# #### Add eqtl results")

# get eQTL's
eqtls = pd.read_table(eqtl_fn)
eqtls.columns = ['eqtl_gname', 'nvar', 'shape1', 'shape2', 'dummy',
                 'sid', 'dist', 'npval', 'slope', 'ppval', 'bpval', 'qval']
print('There are {} eQTLs.'.format(eqtls.shape[0]))

# need to use outer or else you exclude some eQTL's
master = master.merge(eqtls, left_on=['sid', 'gname'], right_on=['sid', 'eqtl_gname'], how='outer')

# add column to filter on eqtl snp status
master['is_eqtl_pair'] = (~master['ppval'].isna()).astype(int)

# add gene names to entries with a missing name (after adding eQTL info)
master.loc[master.gname.isna(), 'gname'] = master.loc[master.gname.isna(), 'eqtl_gname']

# add missing chrA, chrB, startA and startB data for the eQTL rows
master.loc[master.chrA.isna(), 'chrA'] = 'chr' + master.loc[master.chrA.isna(), 'sid'].str.replace(':[0-9]+', '')
master.loc[master.chrB.isna(), 'chrB'] = 'chr' + master.loc[master.chrB.isna(), 'sid'].str.replace(':[0-9]+', '')
master.loc[master.startA.isna(), 'startA'] = (master.loc[master.startA.isna(), 'sid'].str.replace('[0-9]+:', '')).astype(int)
master.loc[master.startA.isna(), 'startA'] -= 1     
master.loc[master.endA.isna(), 'endA'] = master.loc[master.endA.isna(), 'sid'].str.replace('[0-9]+:', '')

print('After outer merging with eqtls master has {} snp-gene pairs.'.format(master.shape[0]))

# #### Add gene meta data 
print("# #### Add gene meta data")
# genes with index as chrom and genename 
query_genes = genes_df.sort_values(['chrom', 'gname']).set_index(['chrom', 'gname'])

def get_gene_meta_from_chrom_gname(query_genes, df, col_idxs=None):
    # add gene positions (for missing gene meta data mostly)
    gene_positions = []
    
    if col_idxs == None:
        for i, sr in df.iterrows():
            gene_info = query_genes.loc[(sr.chrom, sr.gene_name)]

            if len(gene_info) == 0:
                print('Houston, where is my coffee?')
                break
            elif len(gene_info) > 1:
                print('Houston, we have a problem.')
                break
            else:
                gene_positions.append(gene_info.values.tolist()[0])
    else:
        for i, sr in df.iterrows():
            gene_info = query_genes.loc[(sr[col_idxs[0]], sr[col_idxs[1]])]

            if len(gene_info) == 0:
                print(gene_info)
                raise Exception('Houston, where is my coffee?')
                
            elif len(gene_info) > 1:
                #print('Picked the closest gene to the current SNP.')
                dists = np.abs(gene_info['start'].values - sr['startA'])
                closest_idx = np.argmin(dists)
                gene_positions.append(gene_info.values.tolist()[closest_idx][2])                
            else:
                gene_positions.append(gene_info.values.tolist()[0][2])
    return(gene_positions)

gene_ids = get_gene_meta_from_chrom_gname(query_genes, master, col_idxs=[0, 11])
master.loc[:, 'gid'] = gene_ids

# add back the original gene start and end
master = master.merge(orig_genes_df[['start', 'end', 'gene_id', 'strand']], left_on='gid', right_on='gene_id')

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

# #### Add info about closests gene
print("# #### Add info about closests gene")

# check for the closets gene
closets_check = [0] * master.shape[0]
for i, sr in master.iterrows():

    # check closest gene
    rs_gene = (sr.sid, sr.gname)
    if rs_gene in closest_gene.index:
        closets_check[i] = 1

master['is_closest_gene'] = closets_check

# #### Add colocalization data
print("# #### Add colocalization data")

# add colocalization data for SNP and is_coloc_snp columns
tmp_coloc = coloc_sig_full[[
 'pp_H0_Coloc_Summary',
 'pp_H1_Coloc_Summary',
 'pp_H2_Coloc_Summary',
 'pp_H3_Coloc_Summary',
 'pp_H4_Coloc_Summary',
 'rs_id',
 'geneName',
 'ref',
 'alt',
 'AC',
 'AF',
 'AN',
 'slope_gwas',
 'slope_se_gwas',
 'pval_nominal',
 'SampleSize']]
tmp_coloc.rename(columns={'slope_gwas': 'gwas_slope',
                          'slope_se_gwas': 'gwas_slope_se',
                          'pval_nominal': 'gwas_pval_nominal',
                          'geneName': 'gname'}, inplace=True)
master = master.merge(tmp_coloc, on=['rs_id', 'gname'], how='left')

# add column to filter on coloc snp status
master['is_coloc_pair'] = (~master['pp_H4_Coloc_Summary'].isna()).astype(int)

print('After left merging master with the colocalization table there are {} entries.'.format(master.shape[0]))

# #### Add loop data
print('# #### Add loop data')

# check for the loop gene
loop_check = [0] * master.shape[0]
for i, sr in master.iterrows():

    # check closest gene
    rs_gene = (sr.sid, sr.gname)
    if rs_gene in fivekb_gloops_uniq:
        loop_check[i] = 1       

master['has_fithichip_loop'] = loop_check


print('There are {} SNP-Gene loops.'.format(sum(loop_check)))


# #### Add Promoter Capture Chicago Scores

# Include PC HiC data when it's available
if celltype in pchic.columns:
    sg_pchic_view = ['sg_chrA', 'sg_startA', 'sg_endA', 'sg_chrB', 'sg_startB', 'sg_endB', 'chicago_score']
    master = master.merge(sg_pchic_overlaps[sg_pchic_view],
                     left_on=['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB'],
                     right_on=['sg_chrA', 'sg_startA', 'sg_endA', 'sg_chrB', 'sg_startB', 'sg_endB'], how='left')
    master.drop(['sg_chrA', 'sg_startA', 'sg_endA', 'sg_chrB', 'sg_startB', 'sg_endB'], inplace=True, axis=1)


# #### Do the final reordering and saving

# Include PC HiC data when it's available
final_cols = ['sid',
 'rs_id',
 'gname',
 'gid',
 'chrA',
 'endA',    
 'startB',
 'endB',
 'is_eqtl_pair',
 'is_coloc_pair',
 'is_closest_gene',
 'has_fithichip_loop',
 'nvar',
 'shape1',
 'shape2',
 'dist',
 'npval',
 'slope',
 'ppval',
 'bpval',
 'qval',
 'pp_H0_Coloc_Summary',
 'pp_H1_Coloc_Summary',
 'pp_H2_Coloc_Summary',
 'pp_H3_Coloc_Summary',
 'pp_H4_Coloc_Summary',
 'gene_start', 
 'gene_end', 
 'gene_strand', 
 'ref',
 'alt',
 'AC',
 'AF',
 'AN',
 'gwas_slope',
 'gwas_slope_se',
 'gwas_pval_nominal',
 'SampleSize']

if celltype in pchic.columns:
    final_cols += ['chicago_score']
master = master[final_cols]
master.rename(columns={'chrA':'chrom', 'endA': 'snp_pos', 
                       'startB': 'tss_start', 'endB': 'tss_end',
                       'gname': 'gene_name', 'gid': 'gene_id'}, inplace=True)

master.sort_values(['chrom', 'snp_pos', 'tss_start', 'rs_id'], inplace=True)
master.snp_pos = master.snp_pos.astype(int)

# write out the master data
fn = os.path.join(outdir, 'master.tsv')
master.to_csv(fn, sep='\t', header=True, index=False)

fn = os.path.join(outdir, 'master.xlsx')
excel_master = master.sort_values('rs_id').set_index('rs_id')
excel_master.to_excel(fn, na_rep='nan')

#########################################################################################
# Make WashU files 
#########################################################################################

def bedpe_to_WashU_longrange(fn, df):
    """
        Convert from a loop bedpe file into WashU longrange, 
        includes bgzip and tabix of the fn. 
        
        Params
        -------
        fn: str
            path to the longrange output file (without gz)
            
        df: dataframe
            columns 1-6 are as expected and column 7 is the p or q-value. 
            
        Output
        ------
        gzfn: str
            path to the longrange with bgzip compression
        tabix_fn: str
            path to the index of the longrange file
            
    """

    # parsing the data into WashU longrage format
    data = []
    for sr in df.values.tolist():

        # calculate the -log(FDR)
        qval = -np.log(sr[6])

        # get the first pair data
        second_pair_str = '{}:{}-{},{:.5f}'.format(*sr[3:6], qval)
        first_row = sr[0:3] + [second_pair_str]

        # get the second pair data
        first_pair_str = '{}:{}-{},{:.5f}'.format(*sr[0:3], qval)
        second_row = sr[3:6] + [first_pair_str]

        # add each data row
        data.append(first_row)
        data.append(second_row)

    data = sorted(data, key=lambda x: (x[0], x[1], x[2]))

    # writing out the data
    with open(fn, 'w') as f:
        for line in data:
            info = [str(x) for x in line]
            info = '\t'.join(info)
            f.write(info + '\n')
            
    # run bgzip
    cmd = '{} {}'.format(bgzip, fn)
    print(cmd)
    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)

    out, err = job.communicate()
    print('out:', out.decode())
    print('err:', err.decode())
    
    # run tabix
    lrange_gzfn = fn + '.gz'
    cmd = '{} -f {}'.format(tabix, lrange_gzfn)
    print(cmd)
    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)

    out, err = job.communicate()
    print('out:', out.decode())
    print('err:', err.decode())

    print('Created the gzfn: {}'.format(fn + '.gz'))
    print('Created the tabix: {}'.format(fn + '.gz.tbi'))

def bed_WashU_bedgz(fn, df):
    """
        Convert from a bed dataframe into WashU longrange file 
        includes bgzip and tabix of the fn. 
        
        Params
        -------
        fn: str
            path to the longrange output file (without gz)
            
        df: dataframe
            columns 1-3 are as expected and column 7 is the p or q-value. 
            
        Output
        ------
        gzfn: str
            path to the longrange with bgzip compression
        tabix_fn: str
            path to the index of the longrange file
            
    """

    # parsing the data into WashU longrage format
    data = []
    for sr in df.values.tolist():
        data.append(sr[0:4])
    data = sorted(data, key=lambda x: (x[0], x[1], x[2]))

    # writing out the data
    with open(fn, 'w') as f:
        for line in data:
            info = [str(x) for x in line]
            info = '\t'.join(info)
            f.write(info + '\n')
            
    # run bgzip
    cmd = '{} {}'.format(bgzip, fn)
    print(cmd)
    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)

    out, err = job.communicate()
    print('out:', out.decode())
    print('err:', err.decode())
    
    # run tabix
    gzfn = fn + '.gz'
    cmd = '{} -f {}'.format(tabix, gzfn)
    print(cmd)
    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)

    out, err = job.communicate()
    print('out:', out.decode())
    print('err:', err.decode())

    print('Created the gzfn: {}'.format(fn + '.gz'))
    print('Created the tabix: {}'.format(fn + '.gz.tbi'))

def bed_to_WashU_refbed(fn, df):
    """ 
        Convert from a bed dataframe into WashU longrange file 
        includes bgzip and tabix of the fn. 
        
        Params
        -------
        fn: str
            path to the longrange output file (without gz)
            
        df: dataframe
            columns 1-3 are as expected and column 7 is the p or q-value. 
            
        Output
        ------
        gzfn: str
            path to the longrange with bgzip compression
        tabix_fn: str
            path to the index of the longrange file
            
    """

    # parsing the data into WashU longrage format
    data = df.values.tolist()
    data = sorted(data, key=lambda x: (x[0], x[1], x[2]))

    # writing out the data
    with open(fn, 'w') as f:
        for line in data:
            info = [str(x) for x in line]
            info = '\t'.join(info)
            f.write(info + '\n')
    
    # run bgzip
    cmd = '{} -f {}'.format(bgzip, fn) 
    print(cmd)
    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)

    out, err = job.communicate()
    print('out:', out.decode())
    print('err:', err.decode())

    # run tabix
    gzfn = fn + '.gz'
    cmd = '{} {}'.format(tabix, gzfn)
    print(cmd)
    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)

    out, err = job.communicate()
    print('out:', out.decode())
    print('err:', err.decode())

    print('Created the gzfn: {}'.format(fn + '.gz'))
    print('Created the tabix: {}'.format(fn + '.gz.tbi'))


# make the refbed link for genes (status: running)
final_sg_cols = ['chrom', 'gene_start' ,'gene_end', 'gene_name', 'gene_strand']
final_sg_genes = master.loc[(master.has_fithichip_loop == 1), final_sg_cols]

final_sg_genes.gene_start = final_sg_genes.gene_start.astype(int)
final_sg_genes.gene_end = final_sg_genes.gene_end.astype(int)

final_sg_genes = final_sg_genes.loc[~final_sg_genes.duplicated()]
final_sg_genes['chr'] = final_sg_genes['chrom'] 
final_sg_genes['transcript_start'] = final_sg_genes['gene_start']
final_sg_genes['transcript_stop'] = final_sg_genes['gene_end']
final_sg_genes['translation_start'] = final_sg_genes['gene_start']
final_sg_genes['translation_stop'] = final_sg_genes['gene_end']
final_sg_genes['strand'] = final_sg_genes['gene_strand']
final_sg_genes['gene_name'] = final_sg_genes['gene_name']
final_sg_genes['transcript_id'] = final_sg_genes['gene_name']
final_sg_genes['type'] = 'coding'
final_sg_genes['exon_gene_start'] = final_sg_genes['gene_start']
final_sg_genes['exon_stops'] = final_sg_genes['gene_end']
refcols = ['chr', 'transcript_start', 'transcript_stop', 'translation_start',
           'translation_stop', 'strand', 'gene_name', 'transcript_id',
           'type', 'exon_gene_start', 'exon_stops']
final_sg_genes = final_sg_genes.loc[:, refcols]
sg_genes_fn = os.path.join(outdir, 'gs_genes.bed')
bed_to_WashU_refbed(sg_genes_fn, final_sg_genes)


# make the longrange link for raw fithichip data (status: running)
loop_gz = os.path.abspath(loop_fn.replace('.bed', '_WashU.bed.gz'))
loop_tbi = os.path.abspath(loop_fn.replace('.bed', '_WashU.bed.gz.tbi'))
loop_gz_link = os.path.join(outdir, os.path.basename(loop_gz))
loop_tbi_link = os.path.join(outdir, os.path.basename(loop_tbi))
if not os.path.exists(loop_gz_link):
    os.link(loop_gz, loop_gz_link)
    os.link(loop_tbi, loop_tbi_link)


# ## All 5kb Washu Files

# make the longrange link for snp-gene fivekb pairs (status: running)
#fivekb_lrange = fivekb_lrange[fivekb_lrange.rs_id.notna()].reset_index(drop=True)
fivekb_lrange = fivekb_genes.copy()

# convert full for viz
fivekb_lrange = fivekb_lrange[['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB']]
fivekb_lrange.iloc[:, 1] -= 1 
fivekb_lrange['score'] = 0.01

fivekb_snp_gene_pairs_fn = os.path.join(outdir, '5kb.snp_gene_pairs.bed')
bedpe_to_WashU_longrange(fivekb_snp_gene_pairs_fn, fivekb_lrange)

# plot all snp-gene pairs with a loop
fivekb_gloops_lrange = fivekb_gloops.iloc[:, [11,12,13,14,15,16,18]]

fivekb_snp_gene_loops_fn = os.path.join(outdir, '5kb.snp_gene_loops.bed')
bedpe_to_WashU_longrange(fivekb_snp_gene_loops_fn, fivekb_gloops_lrange)


# ## eQTL WashU Files

# make the longrange link for snp-gene eQTL pairs (status: running)
eqtl_lrange = master.loc[master.is_eqtl_pair == 1]
eqtl_lrange = eqtl_lrange[eqtl_lrange.rs_id.notna()].reset_index(drop=True)

# convert full for viz
eqtl_snp_gene_pairs = eqtl_lrange[['chrom', 'snp_pos', 'snp_pos', 'chrom', 'tss_start', 'tss_end', 'gene_name']]
eqtl_snp_gene_pairs.iloc[:, 1] -= 1 
eqtl_snp_gene_pairs['score'] = 0.01

eqtl_snp_gene_pairs_fn = os.path.join(outdir, 'eqtl.snp_gene_pairs.bed')
bedpe_to_WashU_longrange(eqtl_snp_gene_pairs_fn, eqtl_snp_gene_pairs.iloc[:, [0,1,2,3,4,5,7]])

# make the bed for egenes only (status: running)
eqtl_genes = eqtl_lrange.loc[:, ['chrom', 'gene_start', 'gene_end', 'gene_name']]
eqtl_genes_only_fn = os.path.join(outdir, 'eqtl.genes_only.bed')
bed_WashU_bedgz(eqtl_genes_only_fn, eqtl_genes)


# In[49]:


# In[45]:


# make the bed for eSNPs only (status: running)
eqtl_snps = eqtl_lrange.loc[:, ['chrom', 'snp_pos', 'snp_pos', 'rs_id']]
eqtl_snps.iloc[:, 1] -= 1
eqtl_snps = pbt.BedTool.from_dataframe(eqtl_snps)
eqtl_snps = eqtl_snps.slop(b=500, g=gs_fn).to_dataframe()
eqtl_snps_only_fn = os.path.join(outdir, 'eqtl.snps_only.bed')
bed_WashU_bedgz(eqtl_snps_only_fn, eqtl_snps)


# In[50]:


# In[46]:


# # make the longrange link for snp-gene eQTL loops (status: didn't add loop coordinates to master.)
# eqtl_loops = master.loc[master.is_eqtl_pair == 1]


# ## Coloc WashU Files

# In[51]:


# In[47]:


# make the bed for colocalized SNPs (status: running)
final_snps = coloc_sig_df.copy()
final_snps = final_snps.loc[~final_snps.duplicated()]
final_snps = pbt.BedTool.from_dataframe(final_snps)
final_snps = final_snps.slop(b=500, g=gs_fn)
final_snps = final_snps.to_dataframe()[0:4]


# In[52]:


# In[48]:


coloc_snps_only_fn = os.path.join(outdir, 'coloc.snps_only.bed')
bed_WashU_bedgz(coloc_snps_only_fn, final_snps)


# In[53]:


# In[49]:


# make the bed for colocalized genes (status: running)
final_coloc_snp_genes_pairs = master[(master.is_coloc_pair == 1)]
final_coloc_genes = final_coloc_snp_genes_pairs[['chrom', 'gene_start' ,'gene_end', 'gene_name']]
coloc_genes_only_fn = os.path.join(outdir, 'coloc.genes_only.bed')
bed_WashU_bedgz(coloc_genes_only_fn, final_coloc_genes)


# In[54]:


# In[50]:


# make the longrange for colocalized snp-gene pairs (status: running)
final_coloc_snp_genes_pairs_out = final_coloc_snp_genes_pairs[['chrom', 'snp_pos', 'snp_pos', 'chrom', 
                                                               'tss_start', 'tss_end', 'gene_name']]
final_coloc_snp_genes_pairs_out.iloc[:, 1] -= 1
final_coloc_snp_genes_pairs_out.iloc[:, 6] = 0.01

coloc_snp_gene_pairs_fn = os.path.join(outdir, 'coloc.snp_gene_pairs.bed')
bedpe_to_WashU_longrange(coloc_snp_gene_pairs_fn, final_coloc_snp_genes_pairs_out)


# In[55]:


# In[51]:


# make the longrange for colocalized snp-gene loops (status: running)
final_coloc_snp_genes_loops = fivekb_gloops[(fivekb_gloops[10].isin(final_coloc_snp_genes_pairs['rs_id'])) & 
              (fivekb_gloops[12].isin(final_coloc_snp_genes_pairs['gene_id']))]
final_coloc_snp_genes_loops_out = final_coloc_snp_genes_loops.iloc[:, [14,15,16,17,18,19,21]]

coloc_snp_gene_loops_fn = os.path.join(outdir, 'coloc.snp_gene_loops.bed')
bedpe_to_WashU_longrange(coloc_snp_gene_loops_fn, final_coloc_snp_genes_loops_out)


# In[56]:


# In[52]:

final_coloc_snp_genes_loops_out


# In[57]:


# In[53]:


# make the loop anchors as bed files (status: developing)
left = final_coloc_snp_genes_loops[[14,15,16,6]].T.reset_index(drop=True).T
left[3] = 'L-' + left[3] 
right = final_coloc_snp_genes_loops[[17,18,19,6]].T.reset_index(drop=True).T
right[3] = 'R-' + right[3] 
anchors = pd.concat([left, right], ignore_index=True, axis=0)


# In[58]:


# In[54]:
coloc_anchors_fn = os.path.join(outdir, 'coloc.anchors_only.bed')
bed_WashU_bedgz(coloc_anchors_fn, anchors)


# # make the hub json file

#gwas, cline = coloc_fn.split('/')[5:7]

print("# make the hub json file")

gene_refbed_json = {'type': 'refbed',
                    'filename': os.path.basename(sg_genes_fn) + '.gz',
                    'name': 'Gencode V19',
                    'showOnHubLoad': True
                   }

orig_loops_json = {'type': 'longrange',
                   'filename': os.path.basename(loop_gz_link),
                   'name': 'Original Loops',
                   'options': {'displayMode': 'arc', 'color':'red'},
                   'showOnHubLoad': True
                }

# fivekb_snp_gene_pairs_json = {'type': 'longrange',
#                    'filename': os.path.basename(fivekb_snp_gene_pairs_fn) + '.gz',
#                    'name': '5kb SNP-Gene Pairs',
#                    'options': {'displayMode': 'arc', 'color':'purple', 'height': 200},
#                    'showOnHubLoad': False
#                 }

# fivekb_snp_gene_loops_json = {'type': 'longrange',
#                    'filename': os.path.basename(fivekb_snp_gene_loops_fn) + '.gz',
#                    'name': '5kb SNP-Gene Loops',
#                    'options': {'displayMode': 'arc', 'color':'red', 'height': 200},
#                    'showOnHubLoad': False
#                 }


eqtl_snp_gene_pairs_json = {'type': 'longrange',
                 'filename': os.path.basename(eqtl_snp_gene_pairs_fn) + '.gz',
                 'name': 'eQTL SNP-Gene Pairs',
                 'options': {'displayMode': 'arc', 'color':'purple', 'height': 200},
                 'showOnHubLoad': True
                }


eqtl_snps_only_json = {'type': 'bed',
                 'filename': os.path.basename(eqtl_snps_only_fn) + '.gz',
                 'name': 'eQTL SNPs only',
                 'options': {'color':'purple'},
                 'showOnHubLoad': True

                }


eqtl_genes_only_json = {'type': 'bed',
                 'filename': os.path.basename(eqtl_genes_only_fn) + '.gz',
                 'name': 'eQTL genes only',
                 'options': {'color':'purple'},
                 'showOnHubLoad': True
                }


coloc_snp_gene_pairs_json = {'type': 'longrange',
                 'filename': os.path.basename(coloc_snp_gene_pairs_fn) + '.gz',
                 'name': 'coloc SNP-Gene Pairs',
                 'options': {'displayMode': 'arc', 'color':'purple', 'height': 200},
                 'showOnHubLoad': True
                }

coloc_snp_gene_loops_json = {'type': 'longrange',
                 'filename': os.path.basename(coloc_snp_gene_loops_fn) + '.gz',
                 'name': 'coloc SNP-Gene Loops',
                 'options': {'displayMode': 'arc', 'color':'red', 'height': 200},
                 'showOnHubLoad': True
                }

coloc_snps_only_json = {'type': 'bed',
                 'filename': os.path.basename(coloc_snps_only_fn) + '.gz',
                 'name': 'coloc SNPs only',
                 'options': {'color':'purple'},
                 'showOnHubLoad': False
                }


coloc_genes_only_json = {'type': 'bed',
                 'filename': os.path.basename(coloc_genes_only_fn) + '.gz',
                 'name': 'coloc genes only',
                 'options': {'color':'purple'},
                 'showOnHubLoad': False
                }


coloc_anchors_only_json = {'type': 'bed',
                 'filename': os.path.basename(coloc_anchors_fn) + '.gz',
                 'name': 'coloc loop anchors only',
                 'options': {'color':'red'},
                 'showOnHubLoad': False
                }


#hub_json = [orig_loops_json, sg_pairs_json, sg_loops_json, sg_snps_json, sg_genes_json]
hub_json = [gene_refbed_json, 
            #fivekb_snp_gene_pairs_json,
            #fivekb_snp_gene_loops_json,
            eqtl_snps_only_json,
            eqtl_genes_only_json,
            eqtl_snp_gene_pairs_json,
            coloc_snps_only_json,
            coloc_genes_only_json,
            coloc_snp_gene_pairs_json, 
            coloc_anchors_only_json,
            coloc_snp_gene_loops_json,
            orig_loops_json, 
           ]

hub_json_fn = os.path.join(outdir, 'hub.config.json')
with open(hub_json_fn, 'w') as f:
    f.write(json.dumps(hub_json, indent=4))
