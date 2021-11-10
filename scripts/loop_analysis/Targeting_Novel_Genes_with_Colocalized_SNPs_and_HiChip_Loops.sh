#!/usr/bin/env python
# coding: utf-8

# In[1]:


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

# default values for the command line
#sys.argv = [0] * 7
#sys.argv[1] =  'results/main/2021_Nikhil_eQTL/Results/Colocalization/T1D_34012112_Gaulton/'
#sys.argv[1] += 'BLUEPRINT_eQTL_Monocyte/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed'
#sys.argv[2] = 'results/refs/ensembl/gencode.v19.annotation.bed'
#sys.argv[3] = 'results/main/2021_Nikhil_eQTL/Data/FitHiChIP_Loops/CM/FitHiChIP_L/FitHiChIP.interactions_FitHiC_Q0.01.bed'
#sys.argv[4] = 'results/refs/spp/SPP_D-Challenge_networks.xlsx'
#sys.argv[5] = 'results/refs/hg19/hg19.chrom.sizes'
#sys.argv[6] = 'results/main/loop_analysis/washU/'
script="scripts/loop_analysis/Targeting_Novel_Genes_with_Colocalized_SNPs_and_HiChip_Loops.py"


#########################################################################################
# Monocytes #############################################################################
#########################################################################################
echo "# Monocytes"
coloc="results/main/2021_Nikhil_eQTL/Results/Colocalization/T1D_34012112_Gaulton/BLUEPRINT_eQTL_Monocyte/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed"
gencode="results/refs/ensembl/gencode.v19.annotation.bed"
fithichip="results/main/2021_Nikhil_eQTL/Data/FitHiChIP_Loops/CM/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed"
spp="results/refs/spp/SPP_D-Challenge_networks.xlsx"
gs="results/refs/hg19/hg19.chrom.sizes"
outdir="results/main/loop_analysis/CM/"
#python $script $coloc $gencode $fithichip $spp $gs $outdir

#########################################################################################
# Naive B cells #########################################################################
#########################################################################################
echo "# Naive B cells"
coloc="results/main/2021_Nikhil_eQTL/Results/Colocalization/T1D_34012112_Gaulton/DICE_eQTL_B_NAIVE/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed"
gencode="results/refs/ensembl/gencode.v19.annotation.bed"
fithichip="results/main/2021_Nikhil_eQTL/Data/FitHiChIP_Loops/NB/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed"
spp="results/refs/spp/SPP_D-Challenge_networks.xlsx"
gs="results/refs/hg19/hg19.chrom.sizes"
outdir="results/main/loop_analysis/NB/"
#python $script $coloc $gencode $fithichip $spp $gs $outdir

#########################################################################################
# CD4 Naive cells #######################################################################
#########################################################################################
echo "# CD4 Naive cells"
coloc="results/main/2021_Nikhil_eQTL/Results/Colocalization/T1D_34012112_Gaulton/DICE_eQTL_CD4_NAIVE/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed"
gencode="results/refs/ensembl/gencode.v19.annotation.bed"
fithichip="results/main/2021_Nikhil_eQTL/Data/FitHiChIP_Loops/CD4N/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed"
spp="results/refs/spp/SPP_D-Challenge_networks.xlsx"
gs="results/refs/hg19/hg19.chrom.sizes"
outdir="results/main/loop_analysis/CD4N/"
python $script $coloc $gencode $fithichip $spp $gs $outdir

#########################################################################################
# CD8 Naive cells #######################################################################
#########################################################################################
echo "# CD8 Naive cells"
coloc="results/main/2021_Nikhil_eQTL/Results/Colocalization/T1D_34012112_Gaulton/DICE_eQTL_CD8_NAIVE/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed"
gencode="results/refs/ensembl/gencode.v19.annotation.bed"
fithichip="results/main/2021_Nikhil_eQTL/Data/FitHiChIP_Loops/CD8N/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed"
spp="results/refs/spp/SPP_D-Challenge_networks.xlsx"
gs="results/refs/hg19/hg19.chrom.sizes"
outdir="results/main/loop_analysis/CD8N/"
python $script $coloc $gencode $fithichip $spp $gs $outdir

#########################################################################################
# NK CD16POS cells ######################################################################
#########################################################################################
echo "# NK CD16POS cells"
coloc="results/main/2021_Nikhil_eQTL/Results/Colocalization/T1D_34012112_Gaulton/DICE_eQTL_NK_CD16POS/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed"
gencode="results/refs/ensembl/gencode.v19.annotation.bed"
fithichip="results/main/2021_Nikhil_eQTL/Data/FitHiChIP_Loops/NK/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed"
spp="results/refs/spp/SPP_D-Challenge_networks.xlsx"
gs="results/refs/hg19/hg19.chrom.sizes"
outdir="results/main/loop_analysis/NK/"
python $script $coloc $gencode $fithichip $spp $gs $outdir


#########################################################################################
# T follicular helper cells #############################################################
#########################################################################################
echo "# T follicular helper cells"
coloc="results/main/2021_Nikhil_eQTL/Results/Colocalization/T1D_34012112_Gaulton/DICE_eQTL_TFH/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed"
gencode="results/refs/ensembl/gencode.v19.annotation.bed"
fithichip="results/main/2021_Nikhil_eQTL/Data/FitHiChIP_Loops/TFH/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed"
spp="results/refs/spp/SPP_D-Challenge_networks.xlsx"
gs="results/refs/hg19/hg19.chrom.sizes"
outdir="results/main/loop_analysis/TFH/"
python $script $coloc $gencode $fithichip $spp $gs $outdir





#########################################################################################
# Non classical monocytes ###############################################################
#########################################################################################
echo "# T follicular helper cells"
coloc="results/main/2021_Nikhil_eQTL/Results/Colocalization/T1D_34012112_Gaulton/DICE_eQTL_NONCLASSICAL_MONOCYTES/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed"
gencode="results/refs/ensembl/gencode.v19.annotation.bed"
fithichip="results/main/2021_Nikhil_eQTL/Data/FitHiChIP_Loops/NCM/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed"
spp="results/refs/spp/SPP_D-Challenge_networks.xlsx"
gs="results/refs/hg19/hg19.chrom.sizes"
outdir="results/main/loop_analysis/NCM/"
python $script $coloc $gencode $fithichip $spp $gs $outdir
































