#!/bin/bash

#==============================
# Colocalization_Analysis script
#==============================

CodeExec='/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2021_Nikhil_eQTL/Code/Colocalization/Colocalization_Analysis_GWAS_Script.R'



###############################################################
# Setting up eQTL data ########################################
###############################################################

## base directory containing all eQTL data (for easy pathing)
EQTLBASEDIR='/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2021_Nikhil_eQTL/Data/eqtl_sqtl_summ_stats'

## create an array of eQTL samples with a prefix
declare -a eQTLFileList=(${EQTLBASEDIR}"/BLUEPRINT_eQTL/Monocyte.txt.gz" ${EQTLBASEDIR}"/DGN_eQTL/DGN_eQTL.txt" ${EQTLBASEDIR}"/DICE_eQTL/B_NAIVE.txt.gz" ${EQTLBASEDIR}"/GEUVADIS_eQTL/CEU.txt.gz")
declare -a eQTLPrefixList=("BLUEPRINT_eQTL_Monocyte" "DGN_eQTL" "DICE_eQTL_B_NAIVE" "GEUVADIS_eQTL_CEU")

## counting the number of eQTL samples
numeqtlsample=${#eQTLPrefixList[@]}

###############################################################
# Setting up GWAS data ########################################
###############################################################

## base directory containing all T1D GWAS data
GWASDataBaseDir="/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2021_Nikhil_eQTL/Data/T1D_GWAS"

## significant GWAS data (p-value < 5e-8) - file names (only one right now)
GWASFileName='GWAS_input_colocalization_pval_lt_5eMinus8.txt'

## create an array of GWAS samples
declare -a GWASFileList=(${GWASDataBaseDir}"/T1D_25751624/"${GWASFileName})
declare -a GWASPrefixList=("T1D_25751624")

## counting the number of GWAS samples
numgwassample=${#GWASPrefixList[@]}

###############################################################
# Run the colocalization analysis #############################
###############################################################

## BASE Output Directory to contain colocalization results
BASEOUTDIR='/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2021_Nikhil_eQTL/Results/Colocalization'

##=== process individual GWAS files
# for (( i=0; i<${numgwassample}; i++ )); do
for (( i=0; i<1; i++ )); do

    # get the i-th GWAS file
	inp_GWAS_File=${GWASFileList[$i]}
	out_GWAS_Prefix=${GWASPrefixList[$i]}
	echo 'inp_GWAS_File : '${inp_GWAS_File}
	echo 'out_GWAS_Prefix : '${out_GWAS_Prefix}

    # printing the number of GWAS snps
	nline=`cat $inp_GWAS_File | wc -l`
	echo 'number of significant GWAS SNPs for this chromosome : '$nline	

	##=== process individual eQTL files
	# for (( j=0; j<${numeqtlsample}; j++ )); do
	for (( j=0; j<1; j++ )); do

        # get the j-th GWAS file
		inp_eQTL_File=${eQTLFileList[$j]}
		out_eQTL_Prefix=${eQTLPrefixList[$j]}
		echo 'inp_eQTL_File : '${inp_eQTL_File}
		echo 'out_eQTL_Prefix : '${out_eQTL_Prefix}

		##========== output directory for the current execution
		curroutdir=$BASEOUTDIR'/'$out_GWAS_Prefix'/'$out_eQTL_Prefix
		mkdir -p $curroutdir

		scriptfile=`pwd`'/temp_script_Coloc_'${out_GWAS_Prefix}'_'${out_eQTL_Prefix}'.sh'

echo '' > ${scriptfile}
cat <<EOT >> ${scriptfile}
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

Rscript ${CodeExec} ${inp_GWAS_File} ${curroutdir} ${inp_eQTL_File}

EOT

		chmod +x ${scriptfile}
		qsub ${scriptfile}
		
	done 	# end eQTL file list loop
done 	# end GWAS file list loop





