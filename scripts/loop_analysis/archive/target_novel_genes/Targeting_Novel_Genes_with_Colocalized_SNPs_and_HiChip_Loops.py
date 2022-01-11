#!/usr/bin/env python
# coding: utf-8

# Commandline Specification
# Input
# ----------------------------------------------
# 1) colocalization file  x
# 2) gencode gene annotation file in bed format x
# 3) loop data x 
# 4) spp network data x 
# 5) reference chrom sizes x 
# Output
# ----------------------------------------------
# 1) snp-gene loop summary file
# 2) snp-gene pairs longrange file with index
# 3) snp-gene loops longrange file with index

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

# default values for the command line
#sys.argv = [0] * 7
#sys.argv[1] =  'results/main/2021_Nikhil_eQTL/Results/Colocalization/T1D_34012112_Gaulton/'
#sys.argv[1] += 'BLUEPRINT_eQTL_Monocyte/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed'
#sys.argv[2] = 'results/refs/ensembl/gencode.v19.annotation.bed'
#sys.argv[3] = 'results/main/2021_Nikhil_eQTL/Data/FitHiChIP_Loops/CM/FitHiChIP_L/FitHiChIP.interactions_FitHiC_Q0.01.bed'
#sys.argv[4] = 'results/refs/spp/SPP_D-Challenge_networks.xlsx'
#sys.argv[5] = 'results/refs/hg19/hg19.chrom.sizes'
#sys.argv[6] = 'results/main/loop_analysis/washU/'

#########################################################################################
# parsing the commandline arguments #####################################################
#########################################################################################
print("# parsing the commandline arguments #####################################################")
coloc_fn = sys.argv[1]
genes_fn = sys.argv[2]
loop_fn = sys.argv[3]
spp_fn = sys.argv[4]
gs_fn = sys.argv[5]
outdir = sys.argv[6]

# setting the output file names
os.makedirs(outdir, exist_ok=True)
summary_fn = os.path.join(outdir, 'sgloop_summary.xlsx')
sg_pairs_fn = os.path.join(outdir, 'gs_pairs.longrange.bed')
sg_loops_fn = os.path.join(outdir, 'gs_loops.longrange.bed')
sg_genes_fn = os.path.join(outdir, 'gs_genes.refbed.bed')
sg_snps_fn = os.path.join(outdir, 'gs_snps.bed')
hub_json_fn = os.path.join(outdir, 'hub.config.json')


#########################################################################################
# Load the colocalization data ##########################################################
#########################################################################################
print("# Load the colocalization data ##########################################################")

# load the colocalization data
coloc = pd.read_table(coloc_fn)

# extract the most significant according the H4 
coloc_sig_df = coloc[coloc['pp_H4_Coloc_Summary'] > 0.75]
coloc_sig_full = coloc_sig_df.copy(deep=True)
coloc_sig_df = coloc_sig_df[['chr', 'pos', 'rs_id', 'variant_id', 'geneName']]
coloc_sig_df.rename(columns={'pos': 'end'}, inplace=True)
coloc_sig_df.loc[:, 'start'] = coloc_sig_df.loc[:, 'end'] - 1
coloc_sig_df = coloc_sig_df[['chr', 'start', 'end', 'rs_id', 'variant_id', 'geneName']]
coloc_sig_pbt = pbt.BedTool.from_dataframe(coloc_sig_df).sort()

# list of SNP-GENE colocalized pairs 
sg_coloc_set = coloc_sig_full[['rs_id', 'geneName']].values.tolist()
sg_coloc_set = set([tuple(x) for x in sg_coloc_set])
csnp_slop_pbt = coloc_sig_pbt.slop(b=500000, g=gs_fn)

print('In total there are {} colocalized SNPs.'.format(len(coloc_sig_df)))


#########################################################################################
# Load the gene data ####################################################################
#########################################################################################
print("# Load the gene data ####################################################################")

# load the gencode coords
cols = ['chrom', 'start', 'end', 'strand', 'type', 'gene_id', 'gene_name']
gencode = genes_df = pd.read_table(genes_fn, header=None, names=cols)

# extract just the genes
genes_df = gencode.loc[gencode.type.isin(['gene'])]
genes_df.loc[:, 'chrom'] = genes_df['chrom'].astype(str)
genes_df = genes_df.iloc[:, [0,1,2,6,5]]
genes_pbt = pbt.BedTool.from_dataframe(genes_df).sort()


#########################################################################################
# Intersect SNPs and Genes ##############################################################
#########################################################################################
print("# Intersect SNPs and Genes ##############################################################")

# bedtools intersect for SNP (with 500kb slop) and genes
csnps_gene_pairs = csnp_slop_pbt.intersect(genes_pbt, wa=True, wb=True)

# extract the SNP and gene data
csnps_gene_pairs_df = csnps_gene_pairs.to_dataframe()
csnps_gene_pairs_df = csnps_gene_pairs_df.iloc[:, [0,1,2,6,7,8,3,5,9]]

# remove the slop that you added previously. This has to be done before 
# intersecting the snp-gene pairs with the loops
csnps_gene_pairs_df['start'] += 500000
csnps_gene_pairs_df['end'] -= 500000

# converting to bedpe10 plus data format (allows pairtopair to work)
csnps_gene_pairs_df.columns = bedpe_6cols + ['rs_id', 'coloc_gene', 'sg_gene']
csnps_gene_pairs_df['name'] = '.'
csnps_gene_pairs_df['score'] = '.'
csnps_gene_pairs_df['strand1'] = '.'
csnps_gene_pairs_df['strand2'] = '.'
csnps_gene_pairs_df = csnps_gene_pairs_df[bedpe_10cols + ['rs_id', 'coloc_gene', 'sg_gene']]
csnps_gene_pairs_pbt = pbt.BedTool.from_dataframe(csnps_gene_pairs_df)

#########################################################################################
# Loading the Loop Data #################################################################
#########################################################################################
print("# Loading the Loop Data #################################################################")

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

#########################################################################################
# Intersect SNP-Gene Pairs with Loops ###################################################
#########################################################################################
print("# Intersect SNP-Gene Pairs with Loops ###################################################")

# extract loops contain a snp and gene
sgloops_both = csnps_gene_pairs_pbt.pair_to_pair(loops, type='both', **{'is':True})
sgloops_both_df = sgloops_both.to_dataframe(disable_auto_names=True, header=None)

# extract coordinate and meta data columns
cols = list(range(0,6)) + list(range(10,19)) + [20]
sgloops_both_df = sgloops_both_df[cols]

sgloops_both_df.columns = ['chrS', 'startS', 'endS', 'chrG', 'startG', 'endG', 'rs_id', 'coloc_gene', 'sg_gene', 
                           'chrLA', 'startLA', 'endLA', 'chrRA', 'startRA', 'endRA', 'p-value']
                           #'nameL', 'scoreL', 'strandLA', 'strandLB']
sgloops_both_df

## extract loops contain a snp and gene
#sgloops_either = csnps_gene_pairs_pbt.pair_to_pair(loops, type='either', **{'is':True})
#sgloops_either_df = sgloops_either.to_dataframe(disable_auto_names=True, header=None)
#
## extract coordinate and meta data columns
#cols = list(range(0,6)) + list(range(10,19))
#sgloops_either_df = sgloops_either_df[cols]
#sgloops_either_df.columns = ['chrS', 'startS', 'endS', 'chrG', 'startG', 'endG', 'rs_id', 'coloc_gene', 'sg_gene', 
#                      'chrLA', 'startLA', 'endLA', 'chrRA', 'startRA', 'endRA', 'p-value']
#sgloops_either_df

#########################################################################################
# Characterize the GSLoops ##############################################################
#########################################################################################
print("# Characterize the GSLoops ##############################################################")

# loading SPP Differentially expressed genes 
spp = pd.read_excel(spp_fn, sheet_name='3. Experiment DEGs', usecols=[5])
spp_genes = set([x[0] for x in spp.iloc[7:].values])

# getting the SNP-GENE pairs which overlap a fithic loop
sg_hichip_set = sgloops_both_df[['rs_id', 'sg_gene']].values.tolist()
sg_hichip_set = set([tuple(x) for x in sg_hichip_set])

char_data = []
for i, sr in csnps_gene_pairs_df.iterrows():
    
    # check whether this snp-gene pair has a colocalized record
    if (sr.rs_id, sr.sg_gene) in sg_coloc_set:
        coloc_sg = 1
    else:
        coloc_sg = 0
        
    # check whether this snp-gene pair has a sgloop
    if (sr.rs_id, sr.sg_gene) in sg_hichip_set:
        fithichip = 1
    else:
        fithichip = 0
        
    # check whether this gene is differntially expressed according to SPP
    if sr.sg_gene in spp_genes:
        spp = 1
    else:
        spp = 0
    
    char_data.append([coloc_sg, fithichip, spp])

char_data = pd.DataFrame(char_data)
char_data = pd.concat([csnps_gene_pairs_df, char_data], axis=1)
char_data.columns = bedpe_10cols + ['rs_id', 'coloc_egene', 'sg_gene', 'coloc', 'fithichip', 'spp']

final_data = char_data[['rs_id', 'coloc_egene', 'sg_gene', 'coloc', 'fithichip', 'spp']]
final_data = final_data.sort_values(['rs_id', 'coloc', 'fithichip', 'coloc_egene', 'sg_gene'], 
                                    ascending=[True, False, False, True, True])
final_data = final_data.set_index(['rs_id', 'coloc_egene', 'sg_gene'])

# remove duplicate SNP-GENE combination which are tested twice by 
# the colocalization analysis
final_data = final_data[~final_data.reset_index().duplicated(subset=['rs_id', 'sg_gene']).values]
final_data.to_excel(summary_fn, index=True, header=True)

#########################################################################################
# Convert all results to output and WashU visualization results #########################
#########################################################################################
print("# Convert all results to output and WashU visualization results #########################")

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
    print("# parsing the data into WashU longrage format")
    data = []
    if df.shape[1] == 7:
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

    elif df.shape[1] == 6:
        for sr in df.values.tolist():

            # get the first pair data
            second_pair_str = '{}:{}-{},1'.format(*sr[3:6])
            first_row = sr[0:3] + [second_pair_str]

            # get the second pair data
            first_pair_str = '{}:{}-{},1'.format(*sr[0:3])
            second_row = sr[3:6] + [first_pair_str]

            # add each data row
            data.append(first_row)
            data.append(second_row)

    data = sorted(data, key=lambda x: (x[0], x[1], x[2]))

    # writing out the data
    print("# writing out the data")
    with open(fn, 'w') as f:
        for line in data:
            info = [str(x) for x in line]
            info = '\t'.join(info)
            f.write(info + '\n')
            
    # run bgzip
    print("# run bgzip")
    cmd = '{} -f {}'.format(bgzip, fn)
    print(cmd)
    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)

    out, err = job.communicate()
    print('out:', out.decode())
    print('err:', err.decode())
    
    # run tabix
    print("# run tabix")
    lrange_gzfn = fn + '.gz'
    cmd = '{} {}'.format(tabix, lrange_gzfn)
    print(cmd)
    job = sp.Popen(cmd, stderr=sp.PIPE,stdout=sp.PIPE, shell=True)

    out, err = job.communicate()
    print('out:', out.decode())
    print('err:', err.decode())

    print('Created the gzfn: {}'.format(fn + '.gz'))
    print('Created the tabix: {}'.format(fn + '.gz.tbi'))

def bed_to_WashU_bedgz(fn, df):
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


# make the longrange file
print("# make the longrange file for all snp-gene pairs")
final_sg_pairs = csnps_gene_pairs_df.iloc[:, [0,1,2,3,4,5]]
bedpe_to_WashU_longrange(sg_pairs_fn, final_sg_pairs)

# make the longrange file
print("# make the longrange file for all snp-gene loops")
final_sg_loops = sgloops_both_df.iloc[:, [9,10,11,12,13,14,-1]]
bedpe_to_WashU_longrange(sg_loops_fn, final_sg_loops)

# link over the original loops 
print("# link over the original loops")
loop_gz = loop_fn.replace('.bed', '_WashU.bed.gz')
loops_link = os.path.join(outdir, os.path.basename(loop_gz))
if not os.path.exists(loops_link):
    os.link(loop_gz, loops_link)

index_link = os.path.join(outdir, os.path.basename(loop_gz) + '.tbi')
if not os.path.exists(index_link):
    os.link(loop_gz + '.tbi', index_link)

# make the snps bed file 
print("# make the snps bed file")
final_sg_snps = sgloops_both_df.iloc[:, [0,1,2,6]]
final_sg_snps = final_sg_snps.loc[~final_sg_snps.duplicated()]
final_sg_snps = pbt.BedTool.from_dataframe(final_sg_snps)
final_sg_snps = final_sg_snps.slop(b=1000, g=gs_fn)
final_sg_snps = final_sg_snps.to_dataframe()
bed_to_WashU_bedgz(sg_snps_fn, final_sg_snps)

# make the genes refbed file 
print("# make the genes bed file")
final_sg_genes = csnps_gene_pairs_df.iloc[:, [3,4,5,12]]
final_sg_genes = final_sg_genes.loc[~final_sg_genes.duplicated()]
final_sg_genes['strand'] = '+'
final_sg_genes['type'] = 'coding'
sg_genes_cols =  ['chrB', 'startB', 'endB', 'sg_gene', 'strand', 'type']
final_sg_genes = final_sg_genes.loc[:, sg_genes_cols] # done ref code
final_sg_genes = final_sg_genes.iloc[:, [0,1,2,1,2,4,3,3,5,1,2]] # done ref code
bed_to_WashU_refbed(sg_genes_fn, final_sg_genes)

# make the hub json file
print("# make the hub json file")
orig_loops_json = {'type': 'longrange', 
                   'filename': os.path.basename(loops_link),
                   'name': 'Original Loops',
                   'options': {'displayMode': 'arc', 'color':'purple'}
                }

sg_pairs_json = {'type': 'longrange', 
                 'filename': os.path.basename(sg_pairs_fn) + '.gz',
                 'name': 'SNP-Gene Pairs',
                 'options': {'displayMode': 'arc', 'color':'purple'}
                }

sg_loops_json = {'type': 'longrange', 
                 'filename': os.path.basename(sg_loops_fn) + '.gz',
                 'name': 'SNP-Gene Loops', 
                 'options': {'displayMode': 'arc', 'color':'purple'}
                }

sg_snps_json = {'type': 'bed', 
                 'filename': os.path.basename(sg_snps_fn) + '.gz',
                 'name': 'SNP-Gene SNPs',
                 'queryEndpoint': {'name': 'GWAS Catalog', 'endpoint': 'https://www.ebi.ac.uk/gwas/search?query='},
                 'options': {'color':'red'}
                }

sg_genes_json = {'type': 'bed', 
                 'filename': os.path.basename(sg_genes_fn) + '.gz',
                 'name': 'SNP-Gene Genes',
                 'genome': 'Custom',
                 'queryEndpoint': {'name': 'GeneCards', 'endpoint': 'https://www.genecards.org/Search/Keyword?queryString='},
                 'options': {'color':'blue'}
                }

#hub_json = [orig_loops_json, sg_pairs_json, sg_loops_json, sg_snps_json, sg_genes_json]
hub_json = [sg_snps_json, sg_genes_json, sg_loops_json, sg_pairs_json, orig_loops_json]

with open(hub_json_fn, 'w') as f:
    f.write(json.dumps(hub_json, indent=4))


