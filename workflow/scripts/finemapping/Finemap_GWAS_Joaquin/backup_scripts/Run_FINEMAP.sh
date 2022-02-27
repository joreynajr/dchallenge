#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=10GB
#PBS -l walltime=04:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
#source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR


##====================
## script to execute FINEMAP
## using LDStore generated outputs
##====================

##=============
## Parameters - edit
##=============

## total samples (number of columns of genotype)
samplecount=2504	## check the file "all_phase3.psam" in the parent directory (wc -l) to see the number of samples

## finemap executable
finemapexec='/home/sourya/packages/FINEMAP/finemap_v1.4_x86_64/finemap_v1.4_x86_64'

## output directory to store the results for this GWAS data
BaseOutDir_GWAS='/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2021_Nikhil_eQTL/Results/Finemap_GWAS/T1D_32005708'

## this folder will contain the inputs formatted for FINEMAP 
FINEMAPInpDir=$BaseOutDir_GWAS'/FINEMAP_INPUT'

## output directory for FINEMAP when --cond (stepwise conditioning) option is used
BaseOutDir_FINEMAP_cond=$BaseOutDir_GWAS'/FINEMAP_OUTPUT/cond'
mkdir -p $BaseOutDir_FINEMAP_cond

## output directory for FINEMAP when --sss (shotgun stochastic search) option is used
BaseOutDir_FINEMAP_sss=$BaseOutDir_GWAS'/FINEMAP_OUTPUT/sss'
mkdir -p $BaseOutDir_FINEMAP_sss

## number of causal SNPs to be tested
NUMCAUSALSNP=10

## number of threads
NUMTHREAD=8

##=============
## end Parameters - edit
##=============

i=0
while [ 1 == 1 ]; do
	i=`expr $i + 1`
	gwasfile=${BaseOutDir_GWAS}'/GWAS_Regions/Region_'${i}'.txt'
	if [[ ! -f $gwasfile ]]; then
		break
	fi

	## master file containing the finemap execution input
	## when --cond (stepwise conditioning) option is used 
	Masterfile_cond=${BaseOutDir_FINEMAP_cond}'/input_Region_'${i}'.master'
	echo "z;ld;snp;config;cred;log;n_samples" > ${Masterfile_cond}

	## master file containing the finemap execution input
	## when --sss (shotgun stochastic search) option is used 
	Masterfile_sss=${BaseOutDir_FINEMAP_sss}'/input_Region_'${i}'.master'
	echo "z;ld;snp;config;cred;log;n_samples" > ${Masterfile_sss}

	finemap_z_file=${FINEMAPInpDir}'/Region'${i}'.z'
	finemap_ld_file=${FINEMAPInpDir}'/Region'${i}'.ld'

	## construct the master file "Masterfile_cond"
	finemap_out_snp_file=${BaseOutDir_FINEMAP_cond}'/Region'${i}'.snp'
	finemap_out_config_file=${BaseOutDir_FINEMAP_cond}'/Region'${i}'.config'
	finemap_out_cred_file=${BaseOutDir_FINEMAP_cond}'/Region'${i}'.cred'
	finemap_out_log_file=${BaseOutDir_FINEMAP_cond}'/Region'${i}'.log'
	echo ${finemap_z_file}';'${finemap_ld_file}';'${finemap_out_snp_file}';'${finemap_out_config_file}';'${finemap_out_cred_file}';'${finemap_out_log_file}';'${samplecount} >> ${Masterfile_cond}

	## construct the master file "Masterfile_sss"
	finemap_out_snp_file=${BaseOutDir_FINEMAP_sss}'/Region'${i}'.snp'
	finemap_out_config_file=${BaseOutDir_FINEMAP_sss}'/Region'${i}'.config'
	finemap_out_cred_file=${BaseOutDir_FINEMAP_sss}'/Region'${i}'.cred'
	finemap_out_log_file=${BaseOutDir_FINEMAP_sss}'/Region'${i}'.log'
	echo ${finemap_z_file}';'${finemap_ld_file}';'${finemap_out_snp_file}';'${finemap_out_config_file}';'${finemap_out_cred_file}';'${finemap_out_log_file}';'${samplecount} >> ${Masterfile_sss}	

	scriptfile=`pwd`'/temp_FINEMAP_exec_Region_'${i}'.sh'

echo '' > ${scriptfile}
cat <<EOT >> ${scriptfile}
#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=40GB
#PBS -l walltime=10:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
#source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

## fine mapping when --cond (stepwise conditioning) option is used 
${finemapexec} --cond --in-files ${Masterfile_cond} --log --n-causal-snps ${NUMCAUSALSNP}

## fine mapping when --sss (shotgun stochastic search) option is used 
${finemapexec} --sss --in-files ${Masterfile_sss} --log --n-causal-snps ${NUMCAUSALSNP} --n-threads ${NUMTHREAD} 

EOT

	chmod +x ${scriptfile}
	qsub ${scriptfile}

done
