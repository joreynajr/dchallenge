
# Switch to the dchallenge folder
setwd('Y:/Groups/vd-ay/jreyna/projects/dchallenge/')


fn = 'results/main/coloc/Results/Colocalization_SMKN/T1D_34012112_Gaulton/Schmiedel_2018/B-cell_naive/'
fn = paste0(fn, 'chr14_68335157_69335157/merged_SNPs_GWAS_input_colocalization.txt')


merge_SNP_GWAS_Data = read.table(fn, header = T)

##==== quant for dataset1
##====== N: number of SNPs for this gene
dataset1=list(pvalues=merge_SNP_GWAS_Data$pval_nominal,
              type="quant",
              #N=nrow(currloci_GWASdata),
              N=47,
              beta=merge_SNP_GWAS_Data$slope_gwas,
              varbeta=(merge_SNP_GWAS_Data$slope_se_gwas * merge_SNP_GWAS_Data$slope_se_gwas))

## quant for dataset2
## N: number of SNPs for this chromosome (used for colocalization)
dataset2 <- list(pvalues=merge_SNP_GWAS_Data$pvalue,
                 type="quant",
                 #N=nrow(currloci_eqtldata),
                 N=67314,
                 beta=merge_SNP_GWAS_Data$slope_snp,
                 varbeta=(merge_SNP_GWAS_Data$slope_se_snp * merge_SNP_GWAS_Data$slope_se_snp)) #, sdY=stdev_gene_expr)

## colocalization function according to prior probability settings
# using the prior values, p1 and p2 set by the command line
result <- coloc.abf(dataset1=dataset1, dataset2=dataset2, MAF=merge_SNP_GWAS_Data$AF, p12=0.00001)
