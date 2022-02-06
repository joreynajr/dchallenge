#!/bin/bash
set -euo pipefail

#==============================
# Colocalization_Analysis script
#==============================

CodeExec='/mnt/BioAdHoc/Groups/vd-ay/jreyna/projects/dchallenge/results/main/2021_Nikhil_eQTL/Code/Colocalization/Colocalization_Analysis_GWAS_Script_Generalized.R'
setR='/mnt/BioAdHoc/Groups/vd-ay/jreyna/projects/dchallenge/results/main/2021_Nikhil_eQTL/Code/Colocalization/Activate_R-3.6.1.sh'

## base directory containing all eQTL data
echo "## base directory containing all eQTL data"
EQTLBASEDIR='/mnt/BioAdHoc/Groups/vd-ay/jreyna/projects/dchallenge/results/main/2021_Nikhil_eQTL/Data/eqtl_sqtl_summ_stats'

## eQTL samples (complete)
echo "## eQTL samples (complete)"
declare -a eQTLFileList=($(find "${EQTLBASEDIR}/" | grep -v "sQTL" | grep ".gz"))
declare -a eQTLPrefixList=($(find "${EQTLBASEDIR}/" | grep -v "sQTL" | grep ".gz" | awk 'BEGIN{FS="/"; OFS=""} {print $14, "_", gensub(".txt.gz", "", 1, $15)}'))
#for i in ${!eQTLFileList[@]};
#do
#    echo ${eQTLFileList[$i]}	${eQTLPrefixList[$i]}
#done

## number of eQTL samples
numeqtlsample=${#eQTLPrefixList[@]}

## base directory containing all T1D GWAS data
GWASDataBaseDir="/mnt/BioAdHoc/Groups/vd-ay/jreyna/projects/dchallenge/results/main/2021_Nikhil_eQTL/Data/T1D_GWAS"

## significant GWAS data (p-value < 5e-8) - file names
GWASFileName='GWAS_input_colocalization_pval_lt_5eMinus8.txt'

## GWAS samples, Added the Gaulton data
declare -a GWASFileList=("${GWASDataBaseDir}/T1D_34012112_Gaulton/${GWASFileName}")
declare -a GWASPrefixList=("T1D_34012112_Gaulton")

## number of GWAS samples
numgwassample=${#GWASPrefixList[@]}

## BASE Output Directory to contain colocalization results
BASEOUTDIR='/mnt/BioAdHoc/Groups/vd-ay/jreyna/projects/dchallenge/results/main/2021_Nikhil_eQTL/Results/Colocalization_Update'

##=== process individual GWAS files
# for (( i=0; i<${numgwassample}; i++ )); do
for (( i=0; i<1; i++ )); do

	inp_GWAS_File=${GWASFileList[$i]}
	out_GWAS_Prefix=${GWASPrefixList[$i]}
	echo '	inp_GWAS_File : '${inp_GWAS_File}
	echo '	out_GWAS_Prefix : '${out_GWAS_Prefix}

	nline=`cat $inp_GWAS_File | wc -l`
	echo '	number of significant GWAS SNPs for this chromosome : '$nline	

	##=== process individual eQTL files
	# for (( j=0; j<${numeqtlsample}; j++ )); do
	for (( j=0; j<$numeqtlsample; j++ )); do

		inp_eQTL_File=${eQTLFileList[$j]}
		out_eQTL_Prefix=${eQTLPrefixList[$j]}
		echo '		inp_eQTL_File : '${inp_eQTL_File}
		echo '		out_eQTL_Prefix : '${out_eQTL_Prefix}

		##========== output directory for the current execution
		curroutdir=$BASEOUTDIR'/'$out_GWAS_Prefix'/'$out_eQTL_Prefix
		mkdir -p $curroutdir

        scriptdir="$(pwd)/qjobs/"
        mkdir -p $scriptdir
		scriptfile="${scriptdir}/temp_script_Coloc_${out_GWAS_Prefix}_${out_eQTL_Prefix}.sh"

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
cd \$PBS_O_WORKDIR

source $setR
Rscript ${CodeExec} --eqtl-chr {EQTL_CHR} \
                    --eqtl-pos {EQTL_POS} \
                    --eqtl-geneName {EQTL_GENENAME} \
                    --eqtl-dist {EQTL_DIST} \
                    --eqtl-slope {EQTL_SLOPE} \
                    --eqtl-pvalue {EQTL_PVALUE} \
                    --eqtl-FDR {EQTL_FDR} \
                    --eqtl-header {EQTL_HEADER} \
                    ${inp_GWAS_File} \
                    ${curroutdir} \
                    ${inp_eQTL_File}
EOT

		chmod +x ${scriptfile}
		#qsub ${scriptfile}
		
	done 	# end eQTL file list loop
done 	# end GWAS file list loop





