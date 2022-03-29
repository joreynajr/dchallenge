#!/usr/bin/env Rscript
# script to calculate LD SNPs for SNPs in a given SNP file

library(optparse)
library(dplyr)
library(data.table)
library(ggplot2)
options(scipen = 10)
options(datatable.fread.datatable=FALSE)

#########################################################################################
# PLINK parameters - for LD estimation
#########################################################################################
population_list <- c('AFR', 'AMR', 'EAS', 'EUR', 'SAS')
plinkexec <- "/mnt/BioHome/jreyna/software/anaconda3/hichip-db/bin/plink"

# LD threshold of R2
LD_THR <- 0.8

## number of SNPs in LD window
LD_WINDOW_KB <- 1000

## maximum span to compute LD - 1 Mb
LD_WINDOW <- 1000000

#########################################################################################
# Commandline interface 
#########################################################################################
# meant to replace plink_path <- '/mnt/BioAdHoc/Groups/vd-vijay/Ariel/R24_new/LD-snps-1000G/DupsRemoved/'
# meant to replace list_pop_path <- '/mnt/BioAdHoc/Groups/vd-vijay/Ariel/R24_new/LD-snps/lists-pops/'
# meant to replace SNPInfoFile <- \
# paste0('/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2020_IQTL_HiChIP/Data/SNPInfo/SNPInfo_merged_tables/snpinfo_', currchr, '.txt')
option_list = list(
  	make_option(c("--snp-file"), type="character", default=NULL, help="SNP file."),
  	make_option(c("--onekg-dir"), type="character", default=NULL, help="Directory for the 1000G."),
  	make_option(c("--population-dir"), type="character", default=NULL, help="Directory for the 1000G."),
  	make_option(c("--snpinfo-dir"), type="character", default=NULL, help="Directory for the extra SNP information."),
  	make_option(c("--header"), action="store_true", help="Header columns."),
  	make_option(c("--chr-col"), type="integer", default=NULL, help="Chromosome column."),
  	make_option(c("--chr-prefix"), action="store_true", help="Header columns."),
  	make_option(c("--pos-col"), type="integer", default=NULL, help="Position column."),
  	make_option(c("--workdir"), type="character", default=NULL, help="Work directory.")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(sprintf("\n\n**** Commandline Arguments ****"))
print(paste("snp-file:", opt["snp-file"]))
print(paste("onekg-dir:", opt["onekg-dir"]))
print(paste("population-dir:", opt["population-dir"]))
print(paste("snpinfo-dir:", opt["snpinfo-dir"]))
print(paste("header:", opt["header"]))
print(paste("chr-col:", opt["chr-col"]))
print(paste("pos-col:", opt["pos-col"]))
print(paste("workdir:", opt$workdir))
opt["chr-col"] = as.numeric(opt["chr-col"])
opt["pos-col"] = as.numeric(opt["pos-col"])

# make the output directory
dir.create(opt$workdir, showWarnings=F, recursive=T)

#########################################################################################
# Loading the SNP datasets 
#########################################################################################
print("# Loading the SNP datasets")

snp_data <- data.table::fread(as.character(opt["snp-file"]), header=as.logical(opt["header"]))

#colnames(snp_data) <- c("regionID", "GWASLoci", "index", "rsid", "chr",
#                        "pos", "allele1", "allele2", "maf", "beta", "se",
#                        "z", "prob", "log10bf", "mean", "sd", "mean_incl",
#                        "sd_incl","pval")
snp_data <- snp_data %>% rename(chr=as.numeric(opt["chr-col"]),
                                pos=as.numeric(opt["pos-col"]))

if (opt["chr-prefix"] == FALSE){
    snp_data[,"chr"] <- paste0('chr', snp_data[,"chr"])
}

#print(head(snp_data))

cat(sprintf("\n\n=====>> Unique number of snps: %s\n", nrow(snp_data)))

#########################################################################################
# Compute the LD 
#########################################################################################
# This analysis is being done chromosome wise
print("# Compute the LD")

ChrList <- as.vector(sort(unique(snp_data[, 'chr'])))

print('ChrList')
print(ChrList)

cat(sprintf("\n\n====>> List of chromosomes in the snp data: %s ", paste(ChrList, collapse=" ")))

## setting the final output file
Out_Merge_LD_File <- paste0(opt$workdir, '/Out_Merge_LD.txt')
if (file.exists(Out_Merge_LD_File)) {
	system(paste("rm", Out_Merge_LD_File))
}

for (chridx in 1:length(ChrList)) {

    #####################################################################################
	# extract SNPs for the current chromosome 
    #####################################################################################

    # extracting the current chromosome
	currchr <- ChrList[chridx]
	cat(sprintf("\n\nLD for SNPs - considering chromosome : %s ", currchr))

	## setting th 1KG genotype phase information for this chromosome
	bfileinp <- paste0(as.character(opt["onekg-dir"]), currchr, '.1kg.phase3.v5.rsn')

	## SNP information file for this chromosome
	SNPInfoFile <- paste0(opt["snpinfo-dir"], 'snpinfo_', currchr, '.txt')
	SNPInfoData <- data.table::fread(SNPInfoFile, sep=" ", header=T)
	SNPInfoData <- SNPInfoData[, c(1,2,4)]
	colnames(SNPInfoData) <- c('chr', 'pos', 'rsID')

	## snp data for this chromosome
	snp_data_currchr <- snp_data[which(snp_data[,'chr'] == currchr), ]
	print(paste("## snp data for this chromosome - nrow:", nrow(snp_data_currchr)))

	## get the rsID for SNPs from the SNP metadata file
	print("## get the rsID for SNPs from the SNP metadata file")
	snp_data_currchr_sub <- unique(snp_data_currchr[, c('chr', 'pos')])

	#print(colnames(snp_data_currchr_sub))
    #print(colnames(SNPInfoData))

    print(head(snp_data_currchr_sub, 5))
    print(head(SNPInfoData, 5))

    # add RSID information
	snp_data_currchr_sub <- dplyr::inner_join(snp_data_currchr_sub, SNPInfoData)
	print(paste("# add RSID information - nrow:", nrow(snp_data_currchr_sub)))
    if (nrow(snp_data_currchr_sub) == 0){
        print(sprintf('Started with %s SNPs but after merging with RSID there are 0', nrow(snp_data_currchr)))
        print(sprintf('Cannot proceed with %s', currchr))
        next
    }

    #####################################################################################
	# compute LD SNPs
    #####################################################################################
	print("# compute LD SNPs for the SNPs")

    # writing the SNP file for plink
	list_snp_file <- paste0(opt$workdir, '/temp_list_SNPs.txt')
	write.table(snp_data_currchr_sub, list_snp_file, row.names=F, col.names=F, sep="\t", quote=F, append=F)

    # set output and temporary file paths
	outfileprefix <- paste0(opt$workdir, '/temp_snp')
	out_SNP_from_LD_file <- paste0(opt$workdir, '/LD_out_SNP_list.txt')
	out_SNP_from_LD_file_2 <- paste0(opt$workdir, '/LD_out_SNP_list_2.txt')
	
	## loop over each population to compute population-specific LD values
	for (pop_idx in 1:length(population_list)) {

        ## run PLINK with respect to the current population 
		# apply PLINK to get all the SNPs and their LD with respect to the SNPs in list_snp_file
		# use cis-eQTL window of 1 Mb		
		population <- population_list[pop_idx]
		pop_file <- paste0(opt["population-dir"], 'var-super_pop_pop-', population, '.txt')
        cmd = paste0(plinkexec, " --bfile ", bfileinp,
                                " --keep  ", pop_file,
                                " --r2 --ld-snp-list ", list_snp_file,
                                " --ld-window-kb ", LD_WINDOW_KB,
                                " --ld-window ", (LD_WINDOW - 1),
                                " --ld-window-r2 ", LD_THR,
                                " --out ", outfileprefix,
                                "-pop-", population)
        print(cmd)
		system(cmd)

        #################################################################################
        # dump LD pair information
        #################################################################################
        print("# dump LD pair information")
		out_LD_file <- paste0(outfileprefix, "-pop-", population, ".ld")
		if (file.exists(out_LD_file) == TRUE) {
            cmd <- paste0("awk \'{if (NR>1) {print \"chr\"$1\"\t\"$2\"\t\"$3\"\tchr\"$4\"\t\"$5\"\t\"$6\"\t\"$7}}\' ", out_LD_file, " >> ", out_SNP_from_LD_file_2)
			system(cmd)		
		}
	}	# end population loop

    #####################################################################################
    # consolidate LD information across all populations
    #####################################################################################
    print("# consolidate LD information across all populations")
	if (file.exists(out_SNP_from_LD_file_2) == TRUE) {
		## get the unique pairs of rsID along with the highest LD score
		system(paste0("sort -k3,3 -k6,6 -k7,7gr ", out_SNP_from_LD_file_2, " | awk -F\"\t\" \'!seen[$3,$6]++\' - > ", out_SNP_from_LD_file))
	}


    #####################################################################################
	# process the LD information
    #####################################################################################
    if (file.exists(out_SNP_from_LD_file) == TRUE){

        print(paste0("# process LD information for ", currchr))

        LD_SNP_Pair_DF <- data.table::fread(out_SNP_from_LD_file, header=F)
        colnames(LD_SNP_Pair_DF) <- c('chr', 'pos', 'rsID', 'ld_chr', 'ld_pos', 'ld_rsID', 'LD')

        cat(sprintf("\n\n **** nrow LD_SNP_Pair_DF : %s ", nrow(LD_SNP_Pair_DF)))

        ## merge the original SNP with LD information
        mergeDF <- dplyr::inner_join(snp_data_currchr, LD_SNP_Pair_DF)
        cat(sprintf("\n\n ==>> SNPs in LD: %s %s", currchr, nrow(mergeDF)))


        # append the current results to the larger output file
        if (file.exists(Out_Merge_LD_File) == FALSE) {
            write.table(mergeDF, Out_Merge_LD_File, row.names=F, col.names=T, sep="\t", quote=F, append=F)
        } else {
            write.table(mergeDF, Out_Merge_LD_File, row.names=F, col.names=F, sep="\t", quote=F, append=T)
        }
    }

    #####################################################################################
	# remove the temporary files
    #####################################################################################
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

# make an empty file when no LD information has been found 
# for current SNPs
if (file.exists(Out_Merge_LD_File) == FALSE) {
    file.create(Out_Merge_LD_File)
}

