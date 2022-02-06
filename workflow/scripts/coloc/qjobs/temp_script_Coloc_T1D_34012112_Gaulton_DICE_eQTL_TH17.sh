
#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=20GB
#PBS -l walltime=04:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
#source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

source /mnt/BioAdHoc/Groups/vd-ay/jreyna/projects/dchallenge/results/main/2021_Nikhil_eQTL/Code/Colocalization/Activate_R-3.6.1.sh
Rscript /mnt/BioAdHoc/Groups/vd-ay/jreyna/projects/dchallenge/results/main/2021_Nikhil_eQTL/Code/Colocalization/Colocalization_Analysis_GWAS_Script_Generalized.R /mnt/BioAdHoc/Groups/vd-ay/jreyna/projects/dchallenge/results/main/2021_Nikhil_eQTL/Data/T1D_GWAS/T1D_34012112_Gaulton/GWAS_input_colocalization_pval_lt_5eMinus8.txt /mnt/BioAdHoc/Groups/vd-ay/jreyna/projects/dchallenge/results/main/2021_Nikhil_eQTL/Results/Colocalization_Update/T1D_34012112_Gaulton/DICE_eQTL_TH17 /mnt/BioAdHoc/Groups/vd-ay/jreyna/projects/dchallenge/results/main/2021_Nikhil_eQTL/Data/eqtl_sqtl_summ_stats/DICE_eQTL/TH17.txt.gz

