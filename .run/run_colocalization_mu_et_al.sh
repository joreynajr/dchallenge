#########################################################################################
# Getting all GWAS datasets
#########################################################################################

GWASDataBaseDir="/mnt/BioAdHoc/Groups/vd-ay/jreyna/projects/dchallenge/results/main/2021_Nikhil_eQTL/Data/T1D_GWAS"

## significant GWAS data (p-value < 5e-8) - file names
GWASFileName='GWAS_input_colocalization_pval_lt_5eMinus8.txt'

## GWAS samples, Added the Gaulton data
declare -a GWASFileList=("${GWASDataBaseDir}/T1D_34012112_Gaulton/${GWASFileName}")
declare -a GWASPrefixList=("T1D_34012112_Gaulton")

## number of GWAS samples
numgwassample=${#GWASPrefixList[@]}

#########################################################################################
# Getting all eQTL datasets
#########################################################################################

## base directory containing all eQTL data
echo "## base directory containing all eQTL data"
EQTLBASEDIR='/mnt/BioAdHoc/Groups/vd-ay/jreyna/projects/dchallenge/results/main/2021_Nikhil_eQTL/Data/eqtl_sqtl_summ_stats'

## eQTL samples (complete)
echo "## eQTL samples (complete)"
declare -a eQTLFileList=($(find "${EQTLBASEDIR}/" | grep -v "sQTL" | grep ".gz"))
#declare -a eQTLPrefixList=($(find "${EQTLBASEDIR}/" | grep -v "sQTL" | grep ".gz" | awk 'BEGIN{FS="/"; OFS=""} {print $14, "_", gensub(".txt.gz", "", 1, $15)}'))
declare -a eQTLPrefixList=($(find "${EQTLBASEDIR}/" | grep -v "sQTL" | grep ".gz" | awk 'BEGIN{FS="/"; OFS=""} {print $14, gensub(".txt.gz", "", 1, "")}'))

## number of eQTL samples
numeqtlsample=${#eQTLPrefixList[@]}

##=== process individual GWAS files
# for (( i=0; i<${numgwassample}; i++ )); do
outfiles+=""
for (( i=0; i<1; i++ )); do

	#inp_GWAS_File=${GWASFileList[$i]}
	#echo '	inp_GWAS_File : '${inp_GWAS_File}
	out_GWAS_Prefix=${GWASPrefixList[$i]}
	echo '	out_GWAS_Prefix : '${out_GWAS_Prefix}

	##=== process individual eQTL files
	# for (( j=0; j<${numeqtlsample}; j++ )); do
	for (( j=0; j<$numeqtlsample; j++ )); do

		inp_eQTL_File=${eQTLFileList[$j]}
		echo '		inp_eQTL_File : '${inp_eQTL_File}

		out_eQTL_Prefix=${eQTLPrefixList[$j]}
		echo '		out_eQTL_Prefix : '${out_eQTL_Prefix}

        ge_source=$(echo $inp_eQTL_File | xargs basename | cut -d "." -f 1)
        # Adding other files to run with snakemake
        outfiles+=" results/main/coloc/Results/Colocalization_SMK/$out_GWAS_Prefix/$out_eQTL_Prefix/$ge_source/"

        #break

	# end eQTL file list loop
	done
    
    break

# end GWAS file list loop
done

#snakemake --profile workflow/profiles/pbs-torque/ -n $outfiles
#snakemake --profile workflow/profiles/pbs-torque/ $outfiles
#snakemake --profile workflow/profiles/local/ $@ $outfiles
echo "snakemake --profile workflow/profiles/local/ $@ $outfiles"





