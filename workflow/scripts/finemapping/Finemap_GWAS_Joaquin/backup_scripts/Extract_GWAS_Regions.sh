#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=80GB
#PBS -l walltime=40:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
#source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

##====================
## script to extract GWAS significant regions
##====================

##=============
## Parameters - edit
##=============

## directory of this source code
SRCCODEDir='/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2021_Nikhil_eQTL/Code/FineMap/Finemap_GWAS'

## R code for extracting GWAS regions
GWASRegionCodeExec=$SRCCODEDir'/Extract_GWAS_Regions.R'

## input file containing GWAS summary statistics
InpGWASFile='/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2020_IQTL_HiChIP/Data/GWAS_Immune_Studies/T1D_32005708/GCST010681_buildGRCh37.tsv'

## column containing the chromosome name in the GWAS input file
chrcol=1
## column containing the SNP poition in the GWAS input file
poscol=2
## column containing the p-value in the GWAS input file
pvalcol=8

## output directory to store the results for this GWAS data
BaseOutDir_GWAS='/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2021_Nikhil_eQTL/Results/Finemap_GWAS/T1D_32005708'
mkdir -p ${BaseOutDir_GWAS}

##=============
## end Parameters - edit
##=============

## call the script to identify the GWAS regions
Rscript $GWASRegionCodeExec $InpGWASFile ${BaseOutDir_GWAS} $chrcol $poscol $pvalcol



