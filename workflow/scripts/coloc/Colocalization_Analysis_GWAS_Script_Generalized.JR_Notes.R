#!/usr/bin/env Rscript

#========================
# this script is for colocalization between eQTL and GWAS summary statistics
#========================
library(data.table)
library(stringr)
library(coloc)
library(snpStats)
library(ggplot2)
library(ggplotify)
library(argparse)

options(scipen = 10)
options(datatable.fread.datatable=FALSE) 

#########################################################################################
# Setting up command line interface
#########################################################################################
parser <- ArgumentParser(description="Commandline tool for colocalization analyses.")
parser$add_argument('--eqtl-chr', type='integer', default=1)
parser$add_argument('--eqtl-pos', type='integer', default=2)
parser$add_argument('--eqtl-geneName', type='integer', default=3)
parser$add_argument('--eqtl-dist', type='integer', default=4)
parser$add_argument('--eqtl-slope', type='integer', default=5)
parser$add_argument('--eqtl-pvalue', type='integer', default=6)
parser$add_argument('--eqtl-FDR', type='integer', default=NULL)
parser$add_argument('--eqtl-header', type='logical', default=FALSE)
parser$add_argument('--gwas-window', type='integer', default=500000)
parser$add_argument('--p12', type='numeric', default=0.00001)
parser$add_argument('snp_info_dir', metavar='snp-info-dir', type='character')
parser$add_argument('gwas', type='character')
parser$add_argument('eqtl', type='character')
parser$add_argument('mhc-chr6-lb', metavar='mhc_chr6_lb', type='character')
parser$add_argument('mhc-chr6-ub', metavar='mhc_chr6_ub', type='character')
parser$add_argument('outdir', type='character')
params <- parser$parse_args()

# parsing important data and output paths 
inp_GWAS_file <- params$gwas
BaseOutDir <- params$outdir
RefEQTLFile <- params$eqtl

#########################################################################################
# setting parameter values 
#########################################################################################

##======== pp4 threshold for colocalization
##======== not used here - but still worth to note it
THR_POST_PROB <- 0.75

# creating the output directory 
system(paste("mkdir -p", BaseOutDir))

bool_Coloc_Summary_DF <- FALSE

## defining the output summary file which will contain the colocalized SNP-gene pairs
ColocSNPInfoFile <- paste0(BaseOutDir, '/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed')

## defining the log file 
textfile <- paste0(BaseOutDir, '/Out_Summary.log')
outtext <- paste0("\n *** Parameters for colocalization (GWAS) analysis **** \n input GWAS summary statistics file : ",
                    inp_GWAS_file, "\n reference EQTL file : ",
                    RefEQTLFile, "\n BaseOutDir : ", BaseOutDir)
cat(outtext, file=textfile, append=FALSE, sep="\n")

#########################################################################################
## defining a method to estimate SE
#########################################################################################
#=======================
# this function estimates standard error of SNP using beta, MAF values
# references: 1) https://www.biostars.org/p/276869/ (last paragraph), 
# 2) https://www.nature.com/articles/ng.3538, 3) https://www.biostars.org/p/319584/ (last paragraph)
# here b = z / sqrt(2p(1-p)(n+z^2)); where b = beta (regression coefficient), z = z score, p = MAF n = sample size
# so, z can be estimated as  2p(1-p)(b^2)n / sqrt(1 - 2p(1-p)b^2)
# finally standard error is estimated as beta / z
#=======================
estimate_SE <- function(inpdf, beta_column, MAF_column, size_column) { 	
 	b <- inpdf[, beta_column]
 	p <- inpdf[, MAF_column]
 	n <- inpdf[, size_column]
 	numerator <- (2 * p * (1-p) * (b ^ 2) * n)
 	squared_denominator <- (1 - 2 * p * (1-p) * (b ^ 2))

 	Estimated_Z <- rep(0, nrow(inpdf)) 	
 	zero_denom_idx <- which(squared_denominator == 0)
 	nonzero_denom_idx <- setdiff(seq(1,nrow(inpdf)), zero_denom_idx)
 	Estimated_Z[nonzero_denom_idx] <- numerator[nonzero_denom_idx] / sqrt(abs(squared_denominator[nonzero_denom_idx]))
 	Estimated_Z[zero_denom_idx] <- numerator[zero_denom_idx] 	# avoid zero division 

 	# now estimate the standard error (se) which is basically beta / z
 	zero_Z_idx <- which(Estimated_Z == 0)
 	nonzero_Z_idx <- setdiff(seq(1,length(Estimated_Z)), zero_Z_idx)
 	Estimated_SE <- rep(0, length(Estimated_Z))
 	Estimated_SE[nonzero_Z_idx] <- b[nonzero_Z_idx] / Estimated_Z[nonzero_Z_idx]
 	Estimated_SE[zero_Z_idx] <- b[zero_Z_idx]

 	# standard error
 	# outdf <- data.frame(SE=Estimated_SE)
 	# testing with variance 
 	# check https://sciencing.com/calculate-variance-standard-error-6372721.html 	
 	outdf <- data.frame(SE=(Estimated_SE * Estimated_SE * n))

 	return(outdf)
}

#########################################################################################
# compiling the eqtl columns
#########################################################################################
if (is.null(params$eqtl_fdr)){
    eqtl_cols = c(params$eqtl_chr, params$eqtl_pos, params$eqtl_geneName,
                  params$eqtl_dist, params$eqtl_slope, params$eqtl_pvalue)
} else{
    eqtl_cols = c(params$eqtl_chr, params$eqtl_pos, params$eqtl_geneName,
                  params$eqtl_dist, params$eqtl_slope, params$eqtl_pvalue, params$eqtl_fdr)
}

#########################################################################################
## load the GWAS file
#########################################################################################
GWAS_Data <- data.table::fread(inp_GWAS_file, header=T)
if (ncol(GWAS_Data) == 6) {
	colnames(GWAS_Data) <- c('chr', 'pos', 'slope', 'slope_se', 'pval_nominal', 'SampleSize')
} else {
	colnames(GWAS_Data) <- c('chr', 'pos', 'slope', 'slope_se', 'pval_nominal')
}
outtext <- paste0("\n *** Number of reference GWAS SNPs: ", nrow(GWAS_Data))
cat(outtext, file=textfile, append=TRUE, sep="\n")

#########################################################################################
## divide the genome 
#########################################################################################
##=== divide the genome into a fixed set of loci
##=== lead variant: 500 Kb region on both side (1 Mb locus)
##=== split the genome into fixed set of loci

GWASChrList <- as.vector(unique(GWAS_Data[,1]))
for (chridx in 1:length(GWASChrList)) {

    # get the current chromosome
	currchr <- GWASChrList[chridx]	

    #####################################################################################
	## extract the GWAS data for this chromosome
    #####################################################################################

    # get the GWAS summary statistics for the current chromosome
	GWAS_Data_currchr <- GWAS_Data[which(GWAS_Data[,1] == currchr), ]
	outtext <- paste0("\n\n\n *********** \n ==>> processing GWAS chromosome : ",
                      currchr, "\n number of GWAS SNP entries for this chromosome : ",
                      nrow(GWAS_Data_currchr),
                      "\n ****************")
	cat(outtext, file=textfile, append=TRUE, sep="\n")


    #####################################################################################
	## extract the eQTL data for this chromosome
    #####################################################################################

	print("## process the eQTL data for this chromosome")
	## parsing the eQTL data according to commandline specified columns 
	Ref_eQTL_Data <- data.table::fread(RefEQTLFile, select=eqtl_cols, header=params$eqtl_header)
	dump_Ref_eQTL_Data <- Ref_eQTL_Data [, eqtl_cols]

	#print("Printing the Ref_eQTL_Data.")
	#print(head(Ref_eQTL_Data))

    ############### FIX THIS ####################
    #### need to filter eQTL's for data from the current chromosome
	#dump_Ref_eQTL_Data <- Ref_eQTL_Data [, eqtl_cols]
    ############### FIX THIS ####################

    # add FDR data when available
    #print("# adding FDR data when available")
	dump_Ref_eQTL_Data['FDR'] = NA
    if (!is.null(params$eqtl_fdr)){
        dump_Ref_eQTL_Data['FDR'] = Ref_eQTL_Data[, eqtl_fdr] 
    }
    if (is.null(params$eqtl_fdr)){
        dump_Ref_eQTL_Data['FDR'] = NA
    }

	#print("Printing the dump_Ref_eQTL_Data.")
	#print(head(dump_Ref_eQTL_Data))

    # changing the columns names
    print("# adding columns names")
	colnames(dump_Ref_eQTL_Data) <- c('chr', 'pos', 'geneName', 'dist', 'slope', 'pvalue', 'FDR')
	outtext <- paste0("\n ******** Number of reference EQTLs for the current chromosome: ", nrow(dump_Ref_eQTL_Data))
	cat(outtext, file=textfile, append=TRUE, sep="\n")
	if (nrow(dump_Ref_eQTL_Data) == 0) {
		outtext <- paste0("\n !!! As no eQTLs exist, we skip this chromosome !!! ")
		cat(outtext, file=textfile, append=TRUE, sep="\n")
		next
	}

    # setting the filename for the eQTL datset from each dumped chromosome
	temp_eQTL_file <- paste0(BaseOutDir, '/dumped_eQTL_', currchr, '.txt')
	write.table(dump_Ref_eQTL_Data, temp_eQTL_file, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	dump_Ref_eQTL_Data <- data.table::fread(temp_eQTL_file, header=T)		

	########################################################################################
    # add additional SNP information including RS_ID the current chromosome
	########################################################################################
	print("##======== complete SNP information for the current chromosome")

	# merged eQTL with the reference SNP information
    # Adding SNP information including RS ID, ref & alt nucleotide, allele freq and others
	SNPInfoFile <- paste0(params$snp_info_dir, '/snpinfo_', currchr, '.txt')
	SNPInfoData <- data.table::fread(SNPInfoFile, header=T, sep=" ")
	colnames(SNPInfoData) <- c('chr', 'pos', 'variant_id', 'rs_id', 'ref', 'alt', 'AC', 'AF', 'AN')
	cat(sprintf("\n ===>> after reading SNPInfoData "))
	outtext <- paste0("\n ==>>> SNP info file (complete SNP set) : ",
                      SNPInfoFile,
                      "\n ==>> Number of SNPs in the reference SNP information file for the current chromosome : ",
                      nrow(SNPInfoData), " number of columns for this reference dataset : ", ncol(SNPInfoData))
	cat(outtext, file=textfile, append=TRUE, sep="\n")	

	#  merged eQTL with the reference SNP information (by default dplyr will use chr + pos)
	print("##=========== merged eQTL with the reference SNP information")
	merge_EQTL_SNPInfo_DF <- dplyr::inner_join(dump_Ref_eQTL_Data, SNPInfoData)
	cat(sprintf("\n ===>> after creating merge_EQTL_SNPInfo_DF "))
	outtext <- paste0("\n ==>>> number of entries in merge_EQTL_SNPInfo_DF : ", nrow(merge_EQTL_SNPInfo_DF))
	cat(outtext, file=textfile, append=TRUE, sep="\n")

    #####################################################################################
	# extract eQTL fields for colocalization (very important)
    #####################################################################################
	print("##==== extract fields for colocalization")
    cols = c("rs_id", "variant_id", "geneName", "chr", "pos", "dist", "pvalue", "FDR", "slope", "ref", "alt", "AC", "AF", "AN")
	Coloc_CurrChr_DF_SNP <- merge_EQTL_SNPInfo_DF[, cols]

    #####################################################################################
	# estimate standard error using three arguments: beta (column 9), MAF (column 13), size (column 14)
    #####################################################################################
	print("##====== estimate standard error using three arguments: beta (column 9), MAF (column 13), size (column 14)")
	CN <- colnames(Coloc_CurrChr_DF_SNP)
	SE_DF <- estimate_SE(Coloc_CurrChr_DF_SNP, 9, 13, 14)
	Coloc_CurrChr_DF_SNP <- cbind.data.frame(Coloc_CurrChr_DF_SNP, SE_DF$SE)
	colnames(Coloc_CurrChr_DF_SNP) <- c(CN, "slope_se")

    #####################################################################################
	## *** IMPORTANT: excluding SNPs belonging to chr6 MHC locus
    #####################################################################################
    ## within chr6, excluded MHC region from consideration
    ## as suggested in https://www.nature.com/articles/nature22403
	print("##===== *** IMPORTANT: excluding SNPs belonging to chr6 MHC locus")
	mhc_idx <- which((Coloc_CurrChr_DF_SNP[,4] == "chr6") &
                     (Coloc_CurrChr_DF_SNP[,5] >= params$mhc_chr6_lb) &
                     (Coloc_CurrChr_DF_SNP[,5] <= params$mhc_chr6_ub))
	if (length(mhc_idx) > 0) {
		Coloc_CurrChr_DF_SNP <- Coloc_CurrChr_DF_SNP[-mhc_idx, ]
	}
	if (nrow(Coloc_CurrChr_DF_SNP) == 0) {
		next
	}
	outtext <- paste0("\n ==>> Number of SNPs to be used for colocalization for the current",
                      "chromosome (basically number of rows in Coloc_CurrChr_DF_SNP) : ",
                      nrow(Coloc_CurrChr_DF_SNP))
	cat(outtext, file=textfile, append=TRUE, sep="\n")
	
    #####################################################################################
	# define the set of GWAS loci for this chromosome
    #####################################################################################
	print("## define the set of GWAS loci for this chromosome")

	## first, sort the GWAS entries by p-values	
	gwasdata <- GWAS_Data_currchr[order(GWAS_Data_currchr$pval_nominal), ]

	## construct set of GWAS loci for this chromosome
	## after the while loop, "GWAS_Loci_DF" contains the important loci for this chromosome
    ## JR: Note really sure when this loop will terminate
	bool_GWAS_Loci <- FALSE
	while(1) {

		## the top most entry defines the current GWAS loci
		## significant GWAS SNP and 500 Kb in both side
		outtext <- paste0("\n using GWAS SNP - chr : ", currchr, "  pos : ", gwasdata[1, 2])
		cat(outtext, file=textfile, append=TRUE, sep="\n")

        ## 500 Kb on both side of a significant GWAS SNP
		startpos <- max(0, (gwasdata[1, 2] - params$gwas_window))
		endpos <- (gwasdata[1, 2] + params$gwas_window)
		currLociDF <- data.frame(chr=currchr, start=startpos, end=endpos)

        ## initialize the GWAS_Loci_DF or add the new rows
        ## bool_GWAS_Loci helps decide whether to initialize or add new rows
		if (bool_GWAS_Loci == FALSE) {
			GWAS_Loci_DF <- currLociDF
			bool_GWAS_Loci <- TRUE
		} else {
			GWAS_Loci_DF <- rbind.data.frame(GWAS_Loci_DF, currLociDF)
		}

		## discard any GWAS entry belonging to this interval
		idx <- which((gwasdata[,2] >= startpos) & (gwasdata[,2] <= endpos))
		gwasdata <- gwasdata[-idx,]
		if (nrow(gwasdata) == 0) {
			break
		}
	}
	outtext <- paste0("\n ===>> number of rows in GWAS_Loci_DF (GWAS loci for the current chromosome) : ", nrow(GWAS_Loci_DF))	
	cat(outtext, file=textfile, append=TRUE, sep="\n")

    #####################################################################################
	# now process each region/locus with colocalization 
    #####################################################################################
	for (lociidx in 1:nrow(GWAS_Loci_DF)) {

        #################################################################################
        ## extract the eqtl and gwas datasets for the current locus
        #################################################################################

        # extract the start and end positions of the curr region
		startpos <- GWAS_Loci_DF[lociidx, 2]
		endpos <- GWAS_Loci_DF[lociidx, 3]
		outtext <- paste0("\n *** processing GWAS loci : chr : ", currchr, " start : ", startpos, " end : ", endpos)
		cat(outtext, file=textfile, append=TRUE, sep="\n")
		
		## gwas data for the current loci
		currloci_GWASdata <- GWAS_Data_currchr[which((GWAS_Data_currchr[,2] >= startpos) &
                                                     (GWAS_Data_currchr[,2] <= endpos)), ]

		## eQTLs for the current loci
		currloci_eqtldata <- Coloc_CurrChr_DF_SNP[which((Coloc_CurrChr_DF_SNP[,5] >= startpos) &
                                                        (Coloc_CurrChr_DF_SNP[,5] <= endpos)), ]

        # print info about the number of GWAS SNPs
		outtext <- paste0("\n *** Number of GWAS SNPs for this loci : ",
                          nrow(currloci_GWASdata), " number of eQTLs for this loci : ",
                          nrow(currloci_eqtldata))
		cat(outtext, file=textfile, append=TRUE, sep="\n")

        #################################################################################
        ## merge the eqtl and gwas datasets
        #################################################################################
		## merging by chromosome and position
		merge_SNP_GWAS_Data <- merge(currloci_eqtldata, currloci_GWASdata,
                                     by=c("chr", "pos"), all=FALSE, suffixes=c("_snp","_gwas"))
		if (nrow(merge_SNP_GWAS_Data) == 0) {
			outtext <- paste0("\n !!!! merged eQTL and GWAS data has 0 entries -- proceed to the next GWAS loci !!! ")
			cat(outtext, file=textfile, append=TRUE, sep="\n")
			next
		}

        #################################################################################
		# remove NA rows 
        #################################################################################
		NA_idx <- which((is.na(merge_SNP_GWAS_Data$pval_nominal)) |
                        (is.na(merge_SNP_GWAS_Data$pvalue)) |
                        (is.na(merge_SNP_GWAS_Data$slope_snp)) |
                        (is.na(merge_SNP_GWAS_Data$slope_se_snp)) |
                        (is.na(merge_SNP_GWAS_Data$AF)))
		if (length(NA_idx) > 0) {
			merge_SNP_GWAS_Data <- merge_SNP_GWAS_Data[-c(NA_idx), ]
		}

        # skip this loci when no GWAS and eQTL loci are overlapping
		if (nrow(merge_SNP_GWAS_Data) == 0) {
			outtext <- paste0("\n !!!! merged eQTL and GWAS data has 0 entries -- proceed to the next GWAS loci !!! ")
			cat(outtext, file=textfile, append=TRUE, sep="\n")
			next
		}

        # print the number of GWAS-eQTL pairs for this loci
		outtext <- paste0("\n *** Number of entries in merge_SNP_GWAS_Data for this loci : ", nrow(merge_SNP_GWAS_Data))
		cat(outtext, file=textfile, append=TRUE, sep="\n")



        #################################################################################
		# save the input tables for the colocalization analysis 
        #################################################################################

		## dump the input SNPs for colocalization analysis in the specified output directory 
		## for processing the current gene
		CurrLoci_OutDir <- paste0(BaseOutDir, '/', currchr, '_', startpos, '_', endpos)
		system(paste("mkdir -p", CurrLoci_OutDir))

		## write the input GWAS SNP
        snp_gwas_fn = paste0(CurrLoci_OutDir, '/input_GWAS_SNP.txt')
        snp_gwas_df = data.frame(chr=currchr, pos=(GWAS_Loci_DF[lociidx, 2]+GWAS_WINDOW))
        snp_gwas_df = merge(snp_gwas_df, GWAS_Data_currchr)
		write.table(snp_gwas_df, snp_gwas_fn, row.names=F, col.names=T, sep="\t", quote=F, append=F)

		## write the merged GWAS and eQTL data
        gwas_eqtl_fn = paste0(CurrLoci_OutDir, '/merged_SNPs_GWAS_input_colocalization.txt')
		write.table(merge_SNP_GWAS_Data, gwas_eqtl_fn, row.names=F, col.names=T, sep="\t", quote=F, append=F)

        #################################################################################
		# perform the main colocalization routine
        #################################################################################
        # BUG: N is supposed to specify the number of SNPs for this gene but DON'T we need the uniq
        # number of SNPs?

		##==== quant for dataset1
		##====== N: number of SNPs for this gene
		dataset1=list(pvalues=merge_SNP_GWAS_Data$pval_nominal,
                      type="quant",
                      N=nrow(currloci_GWASdata),
                      beta=merge_SNP_GWAS_Data$slope_gwas,
                      varbeta=(merge_SNP_GWAS_Data$slope_se_gwas * merge_SNP_GWAS_Data$slope_se_gwas))

		## quant for dataset2
		## N: number of SNPs for this chromosome (used for colocalization)
		dataset2 <- list(pvalues=merge_SNP_GWAS_Data$pvalue,
                         type="quant",
                         N=nrow(currloci_eqtldata),
                         beta=merge_SNP_GWAS_Data$slope_snp,
                         varbeta=(merge_SNP_GWAS_Data$slope_se_snp * merge_SNP_GWAS_Data$slope_se_snp))	#, sdY=stdev_gene_expr)

		## colocalization function according to prior probability settings
        # using the prior values, p1 and p2 set by the command line
		result <- coloc.abf(dataset1=dataset1, dataset2=dataset2, MAF=merge_SNP_GWAS_Data$AF, p12=params$p12)

        #################################################################################
		## colocalization results: plot sensitivity of the posterior probability computation
        #################################################################################
		if (0) {
			plotfile <- paste0(CurrLoci_OutDir, '/coloc_abf_res_sensitivity_plot.pdf')
			pdf(plotfile, width=8, height=6)
			coloc::sensitivity(result, rule="H4 > 0.5")
			dev.off()
		}

        #################################################################################
        # cleaning up and saving the colocalization results
        #################################################################################

        ## save the summary file
		## result$summary: summary statistic of posterior probability (considering all SNPs)
		coloc_post_prob_summary_file <- paste0(CurrLoci_OutDir, '/coloc_abf_summary_DF1.bed')
		write.table(result$summary, coloc_post_prob_summary_file, row.names=T, col.names=T, sep="\t", quote=F, append=F)

        ## save the results file
		## result$results: detailed statistic of posterior probability per SNP
		## Column 1: SNP.NUM where NUM = index number of SNP - varies from 1 to N
		coloc_post_prob_summary_file_2 <- paste0(CurrLoci_OutDir, '/coloc_abf_results_DF2.bed')
		write.table(result$results, coloc_post_prob_summary_file_2, row.names=F, col.names=T, sep="\t", quote=F, append=F)

		## modify the results (DF2) file to add SNP IDs in addition to the SNP.NUM format
		resDF <- result$results 	
		CN <- colnames(resDF)
		temp_SNP_vec <- as.vector(paste0("SNP.", seq(1,nrow(merge_SNP_GWAS_Data))))
		m <- match(resDF[,1], temp_SNP_vec)
		idx_r <- which(!is.na(m))
		idx_t <- m[!is.na(m)]

		## when data frame merging was done by chr and pos fields
		resDF <- cbind.data.frame(resDF[idx_r, 1],
                                  merge_SNP_GWAS_Data[idx_t, c("rs_id", "pos")],
                                  resDF[idx_r, 2:ncol(resDF)])
		colnames(resDF) <- c(CN[1], 'rs_id', 'pos', CN[2:length(CN)])

        # save the summary file
		coloc_summary_file <- paste0(CurrLoci_OutDir, '/coloc_abf_results_detailed.bed')
		write.table(resDF, coloc_summary_file, row.names=F, col.names=T, sep="\t", quote=F, append=F)
		outtext <- paste0("\n *** creating detailed results -- entries in resDF : ", nrow(resDF), "  matching SNP entries : ", length(idx_r))
		cat(outtext, file=textfile, append=TRUE, sep="\n")

		##======== check if there is a colocalization, w.r.t the PP4 threshold
		##======== in such a case, document the summary statistics

		## extracting the posterior probability (of shared variants) from the summaryy file (CHECK: why not use the result var???)
		pp_H0 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==3) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))
		pp_H1 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==4) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))
		pp_H2 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==5) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))
		pp_H3 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==6) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))
		pp_H4 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==7) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))

        #################################################################################
        # save a final table of colocalized SNPs (PP4 > 0.75)
        #################################################################################
		##====== colocalization - pp.H4 > THR_POST_PROB - potential causal variant
		if ((!is.na(pp_H4)) & (pp_H4 > THR_POST_PROB)) {

            # Extracting SNPs with ppH4 > THR with awk system calls.integer
			coloc_SNP_Set_File <- paste0(CurrLoci_OutDir, '/coloc_SNP_Set.txt')
            cmd <- paste0("awk -F\'[\t]\' \'((NR==1) || (sprintf(\"%0.400f\",$NF)>", THR_POST_PROB, "))\' ", coloc_summary_file, " > ", coloc_SNP_Set_File)
			system(cmd)
			n <- (as.integer(system(paste("cat", coloc_SNP_Set_File, "| wc -l"), intern = TRUE)) - 1)
			outtext <- paste0("\n\n *** COLOCALIZATION -------- number of potential causal variants ( pp.H4 > ", THR_POST_PROB, " : ", n)
			cat(outtext, file=textfile, append=TRUE, sep="\n")

			## extract the SNPs and merge with the eQTL and GWAS statistics 
			tempDF <- data.table::fread(coloc_SNP_Set_File, header=T)
			SNPID <- as.vector(tempDF[, 2])
			SNPPos <- as.vector(tempDF[, 3])
			currDF <- data.frame(chr=rep(currchr, length(SNPID)),
                                 pos=SNPPos,
                                 pp_H0_Coloc_Summary=rep(pp_H0, length(SNPID)),
                                 pp_H1_Coloc_Summary=rep(pp_H1, length(SNPID)),
                                 pp_H2_Coloc_Summary=rep(pp_H2, length(SNPID)),
                                 pp_H3_Coloc_Summary=rep(pp_H3, length(SNPID)),
                                 pp_H4_Coloc_Summary=rep(pp_H4, length(SNPID)))
			mergeDF <- dplyr::inner_join(currDF, merge_SNP_GWAS_Data)

            # initializing the SNPSummaryDF or adding additional rows
			if (bool_Coloc_Summary_DF == FALSE) {						
				SNPSummaryDF <- mergeDF
				bool_Coloc_Summary_DF <- TRUE

			} else {
				SNPSummaryDF <- rbind.data.frame(SNPSummaryDF, mergeDF)
			}
		}

        #################################################################################
		##===== delete temporary objects for each loci
        #################################################################################
		rm(currloci_GWASdata)
		rm(currloci_eqtldata)
		rm(merge_SNP_GWAS_Data)
		rm(dataset1)
		rm(dataset2)
		rm(result)
		rm(resDF)

	}	# end GWAS loci processing loop

    #################################################################################
    ##===== delete temporary objects and files for each chromosome
    #################################################################################
	rm(GWAS_Data_CurrChr)
	rm(dump_Ref_eQTL_Data)
	rm(SNPInfoData)
	rm(merge_EQTL_SNPInfo_DF)
	rm(Coloc_CurrChr_DF_SNP)

	##==== delete temporary files
	if (file.exists(temp_eQTL_file)) {
		system(paste("rm", temp_eQTL_file))
	}

}	# end GWAS chromosome loop

#################################################################################
## save the  the summary statistics
#################################################################################
print("## save the summary statistics")
# dump the final summary file
write.table(SNPSummaryDF, ColocSNPInfoFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)		
