import os 
import sys
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from qmplot import manhattanplot

# collect command line input
gwas_fn = sys.argv[1]
gwas_study = sys.argv[2]
outdir = sys.argv[3]

# test values
#gwas_fn = 'results/main/2021_Nikhil_eQTL/Data/T1D_GWAS/T1D_34012112_Gaulton/GWAS_input_colocalization_pval_lt_5eMinus8.txt'
#gwas_fn = 'results/main/2021_Nikhil_eQTL/Data/T1D_GWAS/T1D_34012112_Gaulton/GWAS_input_colocalization.txt'
#gwas_study = 'T1D_34012112_Gaulton'
#outdir = './'

# make outdir if necessary
os.makedirs(outdir, exist_ok=True)

# load the data
gwas = pd.read_table(gwas_fn)
gwas.columns = ['CHR', 'POS', 'BETA', 'SE', 'P', 'N']
gwas['P'] = gwas['P'].astype(float)

# removing zero value entries 
gwas = gwas.loc[(gwas['P'] != 0),:]

# removing chr in chromosomal name 
gwas['CHR'] = gwas['CHR'].str.replace('chr', '')

# reorganize columns input PLINK format
gwas = gwas[['CHR', 'POS', 'P']]

# graph the manhatton plot
mp_fig = os.path.join(outdir, '{gwas_study}.manhattan_plot.png'.format(gwas_study=gwas_study))
ax = manhattanplot(data=gwas,
                   chrom='CHR',
                   pv='P',
                   xticklabel_kws={"rotation": "vertical"},
                   logp=True, figname=mp_fig)

