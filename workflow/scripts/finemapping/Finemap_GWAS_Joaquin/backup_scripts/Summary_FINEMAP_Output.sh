#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=10GB
#PBS -l walltime=10:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
#source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

CodeExec='/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2021_Nikhil_eQTL/Code/FineMap/Finemap_GWAS/Summary_FINEMAP_Output.R'

## output directory to store the results for this GWAS data
BaseOutDir_GWAS='/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2021_Nikhil_eQTL/Results/Finemap_GWAS/T1D_32005708'

## output directory for FINEMAP when --cond (stepwise conditioning) option is used
BaseOutDir_FINEMAP_cond=$BaseOutDir_GWAS'/FINEMAP_OUTPUT/cond'

## output directory for FINEMAP when --sss (shotgun stochastic search) option is used
BaseOutDir_FINEMAP_sss=$BaseOutDir_GWAS'/FINEMAP_OUTPUT/sss'

## number of causal SNPs to be tested
NUMCAUSALSNP=10

## input GWAS regions
Input_GWAS_Region_File=${BaseOutDir_GWAS}'/GWAS_Regions.txt'

## locuszoom Executable
LocusZoomExec='/mnt/BioAdHoc/Groups/vd-vijay/sourya/Packages/LocusZoom/locuszoom/bin/locuszoom'

## we do not provide LocusZoomExec - it indicates that we won't do locuszoom plots
Rscript $CodeExec ${Input_GWAS_Region_File} ${BaseOutDir_FINEMAP_cond} ${BaseOutDir_FINEMAP_sss} ${NUMCAUSALSNP}

## we provide LocusZoomExec - it indicates that we want to do locuszoom plots
# Rscript $CodeExec ${Input_GWAS_Region_File} ${BaseOutDir_FINEMAP_cond} ${BaseOutDir_FINEMAP_sss} ${NUMCAUSALSNP} ${LocusZoomExec}






