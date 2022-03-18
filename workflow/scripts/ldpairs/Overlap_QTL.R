#!/usr/bin/env Rscript

#======================
# script to et the overlap between IQTLs and the reference set of eQTLs
#======================
library(optparse)
library(dplyr)
library(data.table)
library(ggplot2)
options(scipen = 10)
options(datatable.fread.datatable=FALSE)

##=====
## PLINK parameters - for LD estimation
##=====
population_list <- c('AFR', 'AMR', 'EAS', 'EUR', 'SAS')
plinkexec <- "/home/sourya/packages/PLINK/plink"
plink_path <- '/mnt/BioAdHoc/Groups/vd-vijay/Ariel/R24_new/LD-snps-1000G/'
bcf_path <- paste0(plink_path, 'DupsRemoved/')
list_pop_path <- '/mnt/BioAdHoc/Groups/vd-vijay/Ariel/R24_new/LD-snps/lists-pops/'
LD_THR <- 0.8	# LD threshold of R2

## number of SNPs in LD window
LD_WINDOW_KB <- 1000

## maximum span to compute LD - 1 Mb
LD_WINDOW <- 1000000

#====================================================
option_list = list(
  	make_option(c("--IQTLFile"), type="character", default=NULL, help="IQTL File."),
  	make_option(c("--LoopCol"), type="integer", action="store", default=1, help="Column containing loop information in the IQTL file. Default = 1."),
  	make_option(c("--SNPIDCol"), type="integer", action="store", default=2, help="Column containing SNP ID information in the IQTL file. Default = 2."),
  	make_option(c("--chrCol"), type="integer", action="store", default=3, help="Column containing chromosome information in the IQTL file. Default = 3."),
  	make_option(c("--posCol"), type="integer", action="store", default=4, help="Column containing SNP coordinate in the IQTL file. Default = 4."),
  	make_option(c("--RefAlleleCol"), type="integer", action="store", default=5, help="Column containing information of reference allele in the IQTL file. Default = 5."),
  	make_option(c("--AltAlleleCol"), type="integer", action="store", default=6, help="Column containing information of alternate allele in the IQTL file. Default = 6."),
  	make_option(c("--FDRCol"), type="integer", action="store", default=15, help="Column containing information of FDR (RAQUAL default model) in the IQTL file. Default = 15."),
  	make_option(c("--EQTLFile"), type="character", default=NULL, help="Reference eQTL File."),
  	make_option(c("--EQTLSNPIDCol"), type="integer", action="store", default=4, help="Column containing SNP ID in the reference eQTL file. Default = 4."),
  	make_option(c("--EQTLGeneIDCol"), type="integer", action="store", default=3, help="Column containing gene ID information in the reference eQTL file. Default = 3."),
  	make_option(c("--EQTLchrCol"), type="integer", action="store", default=1, help="Column containing chromosome information in the reference eQTL file. Default = 1."),
  	make_option(c("--EQTLposCol"), type="integer", action="store", default=2, help="Column containing SNP coordinate in the reference eQTL file. Default = 2."),
  	make_option(c("--EQTLGeneNameCol"), type="integer", action="store", default=11, help="Column containing gene name in the reference eQTL file. Default = 11."),
  	make_option(c("--EQTLPValCol"), type="integer", action="store", default=5, help="Column containing p-value in the reference eQTL file. Default = 5."),
  	make_option(c("--EQTLFDRCol"), type="integer", action="store", default=6, help="Column containing FDR in the reference eQTL file. Default = 6."),
  	make_option(c("--OutDir"), type="character", default=NULL, help="Output directory.")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

system(paste("mkdir -p", opt$OutDir))
logfile <- paste0(opt$OutDir, '/out.log')
sink(logfile)

## directory to contain overlapping eQTL
OverlapOutDir <- paste0(opt$OutDir, '/Overlap_eQTL')
system(paste("mkdir -p", OverlapOutDir))

## directory to contain LD of eQTL
LDOutDir <- paste0(opt$OutDir, '/LD_eQTL')
system(paste("mkdir -p", LDOutDir))

cat(sprintf("\n\n **** Input parameters **** \n IQTLFile : %s \n LoopCol : %s \n SNPIDCol : %s \n chrCol : %s \n posCol : %s \n RefAlleleCol : %s \n AltAlleleCol : %s \n FDRCol : %s \n EQTLFile : %s \n EQTLSNPIDCol : %s \n EQTLGeneIDCol : %s \n EQTLchrCol : %s \n EQTLposCol : %s \n EQTLGeneNameCol : %s \n EQTLPValCol : %s \n EQTLFDRCol : %s \n OutDir : %s ", opt$IQTLFile, opt$LoopCol, opt$SNPIDCol, opt$chrCol, opt$posCol, opt$RefAlleleCol, opt$AltAlleleCol, opt$FDRCol, opt$EQTLFile, opt$EQTLSNPIDCol, opt$EQTLGeneIDCol, opt$EQTLchrCol, opt$EQTLposCol, opt$EQTLGeneNameCol, opt$EQTLPValCol, opt$EQTLFDRCol, opt$OutDir))

## read IQTL data
IQTLData <- data.table::fread(opt$IQTLFile, header=T)
IQTLData <- IQTLData[, c(opt$LoopCol, opt$SNPIDCol, opt$chrCol, opt$posCol, opt$RefAlleleCol, opt$AltAlleleCol, opt$FDRCol)]
colnames(IQTLData) <- c('Loop', 'SNPID', 'chr', 'pos', 'Ref', 'Alt', 'FDR_RASQ_Def')

## read reference eQTL data
eQTLData <- data.table::fread(opt$EQTLFile, header=T)
eQTLData <- eQTLData[, c(opt$EQTLchrCol, opt$EQTLposCol, opt$EQTLSNPIDCol, opt$EQTLGeneIDCol, opt$EQTLGeneNameCol, opt$EQTLPValCol, opt$EQTLFDRCol)]
colnames(eQTLData) <- c('chr', 'pos', 'rsID', 'eQTL_geneID', 'eQTL_geneName', 'eQTL_pval', 'eQTL_FDR')

##=====================
## compute the overlapping eQTL and IQTLs
## merge by SNP coordinates
##=====================
mergeDF <- dplyr::inner_join(IQTLData, eQTLData)
write.table(mergeDF, paste0(OverlapOutDir, '/Merged_IQTL_eQTL.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)

## dump unique SNPs
IQTL_SNP_DF <- unique(IQTLData[, c('chr', 'pos', 'SNPID')])
write.table(IQTL_SNP_DF, paste0(OverlapOutDir, '/Input_IQTL_Unique_SNPs.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)
eQTL_SNP_DF <- unique(eQTLData[, c('chr', 'pos', 'rsID')])
write.table(eQTL_SNP_DF, paste0(OverlapOutDir, '/Input_eQTL_Unique_SNPs.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)
merged_IQTL_eQTL_SNP_DF <- unique(mergeDF[, c('chr', 'pos', 'rsID')])
write.table(merged_IQTL_eQTL_SNP_DF, paste0(OverlapOutDir, '/Overlap_IQTL_eQTL_Unique_SNPs.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)

## also write the IQTLs which are not overlapping with eQTLs
write.table(dplyr::anti_join(IQTLData, unique(mergeDF[, c('chr', 'pos')])), paste0(OverlapOutDir, '/IQTL_Not_Overlapping_with_eQTL.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)

## dump output logs
cat(sprintf("\n\n =====>> Overlap of eQTL and IQTL \n\n -->> Number of unique SNPs in IQTL data : %s \n -->> Number of unique SNPs in eQTL data : %s \n -->> Number of unique SNPs in merged eQTL + IQTL data : %s ", nrow(IQTL_SNP_DF), nrow(eQTL_SNP_DF), nrow(merged_IQTL_eQTL_SNP_DF)))

##=====================
## compute the LD between eQTL and IQTLs
## chromosome wise
##=====================
## final output file
Out_Merge_IQTL_eQTL_LD_File <- paste0(LDOutDir, '/Out_Merge_IQTL_eQTL_LD.txt')
if (file.exists(Out_Merge_IQTL_eQTL_LD_File)) {
	system(paste("rm", Out_Merge_IQTL_eQTL_LD_File))
}

## also we store the IQTL SNPs and eQTL SNPs which are in LD with each other
Out_IQTL_LD_File <- paste0(LDOutDir, '/Out_IQTL_LD_with_eQTL.txt')
Out_eQTL_LD_File <- paste0(LDOutDir, '/Out_eQTL_LD_with_IQTL.txt')
Out_IQTL_Not_in_LD_File <- paste0(LDOutDir, '/Out_IQTL_Not_LD_with_eQTL.txt')

ChrList <- as.vector(sort(unique(IQTLData[, c('chr')])))
cat(sprintf("\n\n =====>> List of chromosomes in the IQTL data : %s ", paste(ChrList, collapse=" ")))

for (chridx in 1:length(ChrList)) {
	currchr <- ChrList[chridx]
	cat(sprintf("\n\n LD between eQTLs and IQTLs - considering chromosome : %s ", currchr))
	## 1000G genotype phase information for this chromosome
	bfileinp <- paste0(bcf_path, currchr, '.1kg.phase3.v5.rsn')
	## SNP information file for this chromosome
	SNPInfoFile <- paste0('/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2020_IQTL_HiChIP/Data/SNPInfo/SNPInfo_merged_tables/snpinfo_', currchr, '.txt')
	SNPInfoData <- data.table::fread(SNPInfoFile, header=T)
	SNPInfoData <- SNPInfoData[, c(1,2,4)]
	colnames(SNPInfoData) <- c('chr', 'pos', 'rsID')

	## eQTL data for this chromosome
	eQTLData_currchr <- eQTLData[which(eQTLData[,1] == currchr), ]
	colnames(eQTLData_currchr) <- c('eQTL_chr', 'eQTL_pos', 'eQTL_rsID', 'eQTL_geneID', 'eQTL_geneName', 'eQTL_pval', 'eQTL_FDR')
	## IQTL data for this chromosome
	IQTLData_currchr <- IQTLData[which(IQTLData[,3] == currchr), ]
	colnames(IQTLData_currchr) <- c('Loop', 'SNPID', 'chr', 'pos', 'Ref', 'Alt', 'FDR_RASQ_Def')

	## get the rsID for IQTL SNPs from the SNP information file
	IQTLData_currchr_sub <- unique(IQTLData_currchr[, c(3,4)])
	IQTLData_currchr_sub <- dplyr::inner_join(IQTLData_currchr_sub, SNPInfoData)









	##============
	## compute LD SNPs for the IQTL SNPs
	##============ 
	list_snp_file <- paste0(LDOutDir, '/temp_list_IQTL_SNPs.txt')
	write.table(IQTLData_currchr_sub, list_snp_file, row.names=F, col.names=F, sep="\t", quote=F, append=F)

	## run PLINK
	outfileprefix <- paste0(LDOutDir, '/temp_snp')
	out_SNP_from_LD_file <- paste0(LDOutDir, '/LD_out_SNP_list.txt')
	out_SNP_from_LD_file_2 <- paste0(LDOutDir, '/LD_out_SNP_list_2.txt')
	if (file.exists(out_SNP_from_LD_file)) {
		system(paste("rm", out_SNP_from_LD_file))
	}
	if (file.exists(out_SNP_from_LD_file_2)) {
		system(paste("rm", out_SNP_from_LD_file_2))
	}
	
	## population loop to compute LD
	for (pop_idx in 1:length(population_list)) {
		# with respect to the current population 
		# apply PLINK to get all the SNPs and their LD with respect to the SNPs in list_snp_file
		# use cis-eQTL window of 1 Mb		
		population <- population_list[pop_idx]
		pop_file <- paste0(list_pop_path, 'var-super_pop_pop-', population, '.txt')
        cmd = paste0(plinkexec, " --bfile ", bfileinp,
                                " --keep  ", pop_file,
                                " --r2 --ld-snp-list ", list_snp_file,
                                " --ld-window-kb ", LD_WINDOW_KB,
                                " --ld-window ", (LD_WINDOW - 1),
                                " --ld-window-r2 ", LD_THR,
                                " --out ", outfileprefix,
                                "-pop-", population)
		system(cmd)
        exit()
		## output file containing LD information
		out_LD_file <- paste0(outfileprefix, "-pop-", population, ".ld")
		if (file.exists(out_LD_file) == TRUE) {
			## dump LD pair information
			system(paste0("awk \'{if (NR>1) {print \"chr\"$1\"\t\"$2\"\t\"$3\"\tchr\"$4\"\t\"$5\"\t\"$6\"\t\"$7}}\' ", out_LD_file, " >> ", out_SNP_from_LD_file_2))		
		}
	}	# end population loop
	if (file.exists(out_SNP_from_LD_file_2) == TRUE) {
		## consolidate LD information across all populations
		## get the unique pairs of rsID along with the highest LD score
		system(paste0("sort -k3,3 -k6,6 -k7,7gr ", out_SNP_from_LD_file_2, " | awk -F\"\t\" \'!seen[$3,$6]++\' - > ", out_SNP_from_LD_file))
		# system(paste("rm", out_SNP_from_LD_file_2))
	}


















	## read the LD information
	LD_SNP_Pair_DF <- data.table::fread(out_SNP_from_LD_file, header=F)
	colnames(LD_SNP_Pair_DF) <- c('chr', 'pos', 'rsID', 'eQTL_chr', 'eQTL_pos', 'eQTL_rsID', 'LD')
	cat(sprintf("\n\n **** nrow LD_SNP_Pair_DF : %s ", nrow(LD_SNP_Pair_DF)))

	## merge with both IQTL and eQTL information
	mergeDF <- dplyr::inner_join(dplyr::inner_join(IQTLData_currchr, LD_SNP_Pair_DF), eQTLData_currchr)
	cat(sprintf("\n\n ==>> SNPs in LD between IQTL and eQTL - for chromosome - %s - total number of entries : %s ", currchr, nrow(mergeDF)))

	## get the unique IQTL SNPs
	Unique_IQTL_DF <- unique(mergeDF[, c('Loop', 'SNPID', 'chr', 'pos', 'Ref', 'Alt', 'FDR_RASQ_Def')])

	## get the unique eQTL SNPs
	Unique_eQTL_DF <- unique(mergeDF[, c('eQTL_chr', 'eQTL_pos', 'eQTL_rsID', 'eQTL_geneID', 'eQTL_geneName', 'eQTL_pval', 'eQTL_FDR')])

	## get the IQTLs not in LD with eQTLs
	IQTL_Not_LD_eQTL_DF <- dplyr::anti_join(IQTLData_currchr, unique(Unique_IQTL_DF[, c('chr', 'pos')]))

	if (file.exists(Out_Merge_IQTL_eQTL_LD_File) == FALSE) {
		write.table(mergeDF, Out_Merge_IQTL_eQTL_LD_File, row.names=F, col.names=T, sep="\t", quote=F, append=F)
		write.table(Unique_IQTL_DF, Out_IQTL_LD_File, row.names=F, col.names=T, sep="\t", quote=F, append=F)
		write.table(Unique_eQTL_DF, Out_eQTL_LD_File, row.names=F, col.names=T, sep="\t", quote=F, append=F)
		write.table(IQTL_Not_LD_eQTL_DF, Out_IQTL_Not_in_LD_File, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	} else {
		write.table(mergeDF, Out_Merge_IQTL_eQTL_LD_File, row.names=F, col.names=F, sep="\t", quote=F, append=T)
		write.table(Unique_IQTL_DF, Out_IQTL_LD_File, row.names=F, col.names=F, sep="\t", quote=F, append=T)
		write.table(Unique_eQTL_DF, Out_eQTL_LD_File, row.names=F, col.names=F, sep="\t", quote=F, append=T)		
		write.table(IQTL_Not_LD_eQTL_DF, Out_IQTL_Not_in_LD_File, row.names=F, col.names=F, sep="\t", quote=F, append=T)
	}

	## remove the temporary files
	system(paste0("rm ", outfileprefix, "*"))

	if (file.exists(list_snp_file)) {
		system(paste("rm", list_snp_file))
	}	
	if (file.exists(out_SNP_from_LD_file)) {
		system(paste("rm", out_SNP_from_LD_file))
	}
	if (file.exists(out_SNP_from_LD_file_2)) {
		system(paste("rm", out_SNP_from_LD_file_2))
	}
}



