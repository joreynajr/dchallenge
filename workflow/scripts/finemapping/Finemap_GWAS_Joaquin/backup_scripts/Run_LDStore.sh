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

##====================
## script to execute LDSTORE
## from the input GWAS summary statistics
## and using the PLINK formatted genotype data
##====================

##=============
## Parameters - edit
##=============

## directory of this source code
SRCCODEDir='/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2021_Nikhil_eQTL/Code/FineMap/Finemap_GWAS'

## R code for executing preliminery check before LDStore
LDStoreTestCodeExec=$SRCCODEDir'/Check_SNPs_Input_Ref.R'

## ld store executable
ldstoreexec='/home/sourya/packages/ldstore/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64'

## base directory containing genotyping information
## in the PLINK formatted files, from 1000G SNPs
## to be used for fine mapping GWAS 
## sample count - indicating the sample size

# GENOTYPEDIR='/mnt/BioAdHoc/Groups/vd-vijay/Cristian/sourya/R24_genotyping'
# samplecount=176	## for DICE genotype

# GENOTYPEDIR='/mnt/BioAdHoc/Groups/vd-vijay/Ariel/R24_new/LD-snps-1000G'
GENOTYPEDIR='/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/1000G_Phase3_Genotype/ALL'
samplecount=2504	## check the file "all_phase3.psam" in the parent directory (wc -l) to see the number of samples

## column containing the chromosome name in the GWAS input file
chrcol=1
## column containing the SNP poition in the GWAS input file
poscol=2
## column containing the p-value in the GWAS input file
pvalcol=8
# ## column containing the allele 1 information in the GWAS input file 
# allele1col=3
# ## column containing the other allele (allele 2) information in the GWAS input file
# allele2col=4
## column containing the allele frequency information in the GWAS input file
afcol=5
## column containing the beta information in the GWAS input file
betacol=6
## column containing the standard error information in the GWAS input file
secol=7

## output directory to store the results for this GWAS data
BaseOutDir_GWAS='/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2021_Nikhil_eQTL/Results/Finemap_GWAS/T1D_32005708'
mkdir -p ${BaseOutDir_GWAS}

## this folder will contain the inputs formatted for FINEMAP 
FINEMAPInpDir=$BaseOutDir_GWAS'/FINEMAP_INPUT'
mkdir -p ${FINEMAPInpDir}

## output log file
logfile=${BaseOutDir_GWAS}'/out.log'
echo 'output log' > $logfile

##=============
## end Parameters - edit
##=============

##=======================
## Step 2: execute LDStore
##=======================

## create a temporary file
tempfile=${FINEMAPInpDir}'/t1.txt'

##=== process individual GWAS regions
i=0
while [ 1 == 1 ]; do
	i=`expr $i + 1`
	gwasfile=${BaseOutDir_GWAS}'/GWAS_Regions/Region_'${i}'.txt'
	if [[ ! -f $gwasfile ]]; then
		break
	fi
	echo 'processing gwasfile : '$gwasfile >> $logfile
	finemap_z_file=${FINEMAPInpDir}'/Region'${i}'.z'

	## master file which will be used as the input of FINEMAP
	Masterfile=${FINEMAPInpDir}'/input.Region_'${i}'.master'
	echo "z;bgen;bgi;bcor;ld;n_samples;bdose" > ${Masterfile}

	## get the chromosome
	awk '{if (NR==2) {print $0}}' ${gwasfile} > $tempfile	
	currchr=`awk '{if (NR==1) {print $1}}' ${tempfile}`
	echo 'current chromosome : '$currchr >> $logfile
	## get the chromosome number
	sed -i 's/chr//g' ${tempfile}
	currchrnum=`awk '{if (NR==1) {print $1}}' ${tempfile}`
	echo 'current chromosome number : '$currchrnum >> $logfile

	## reference 1000G genotypes - variants are referred by chrNum:pos
	## we are using all the reference panels
	bgenfile=$GENOTYPEDIR'/all_phase3_1000G_'$currchr'.bgen'
	bgenbgifile=$GENOTYPEDIR'/all_phase3_1000G_'$currchr'.bgen.bgi'
	bimfile=$GENOTYPEDIR'/all_phase3_1000G_'$currchr'.bim'

	## use chromosome:position format for variant ID
	## this file will also be useful for FINEMAP
	## note: 6th field requires MAF - so use < 0.5 and convert the allele frequency if required
	## space separated file
	echo "rsid chromosome position maf beta se" > ${finemap_z_file}
	awk -v n="${currchrnum}" -v C="$chrcol" -v S="$poscol" -v F="$afcol" -v B="$betacol" -v E="$secol" '{if (NR>1) {if ($F>0.5) {f=(1-$F)} else {f=$F}; print n":"$S" "n" "$S" "f" "$B" "$E}}'  $gwasfile >> ${finemap_z_file}	

	##========= before running LDStore, check if these SNPs are present in the reference genotype (PLINK) output
	##========= otherwise filter those entries
	n=`cat ${finemap_z_file} | wc -l`
	echo 'before filtering SNPs from FINEMAP input - number of entries : '$n >> $logfile
	Rscript ${LDStoreTestCodeExec} ${finemap_z_file} ${bimfile}
	n=`cat ${finemap_z_file} | wc -l`
	echo 'after filtering SNPs from FINEMAP input - number of entries : '$n >> $logfile

	##======== another filtering - sourya
	##======== sometimes same SNPs (same rsID, position) exist multiple times, even in the original GWAS data - filter them
	awk -F" " '!seen[$1]++' ${finemap_z_file} > $tempfile
	mv $tempfile ${finemap_z_file}
	n=`cat ${finemap_z_file} | wc -l`
	echo 'after removing any duplicate SNPs from FINEMAP input - number of entries : '$n >> $logfile

	##========= write the master file entries
	outbcorfile=${FINEMAPInpDir}'/Region'${i}'.bcor'
	outldfile=${FINEMAPInpDir}'/Region'${i}'.ld'
	bdosefile=${FINEMAPInpDir}'/Region'${i}'.bdose'
	echo ${finemap_z_file}';'${bgenfile}';'${bgenbgifile}';'${outbcorfile}';'${outldfile}';'${samplecount}';'${bdosefile} >> ${Masterfile}

	scriptfile=`pwd`'/temp_LDStore_exec_Region_'${i}'.sh'

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

## now execute ldstore using the generared master file
${ldstoreexec} --in-files ${Masterfile} --read-only-bgen --write-text --write-bdose --n-threads 8 --memory 40

EOT

	chmod +x ${scriptfile}
	qsub ${scriptfile}

done

