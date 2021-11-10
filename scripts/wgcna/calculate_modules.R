library(WGCNA)
library(sets)
library(biomaRt)
library(fgsea)
library(anRichment)

allowWGCNAThreads()
enableWGCNAThreads()
options(stringsAsFactors = FALSE);
outdir = 'results/main/'

##########################################################
########### 1. Data input and cleaning ###################
##########################################################

# load RNA-seq data (samples by genes)
fn = "results/refs/datasets/GSE123658/GSE123658_read_counts.gene_level.txt"
geneExpr = read.csv(fn, sep="\t");

# tidy up the RNA-seq data into a new variable 
rownames(geneExpr) = geneExpr$Samples; #set the row names
geneExpr = as.data.frame(geneExpr[, -c(1)]); # remove the first column
geneExpr = t(geneExpr) 

# load gene length data (genes + lengths)
gene_info = read.table('results/refs/ensembl/Homo_sapiens.GRCh38.104.bed')
colnames(gene_info) = c('chrom', 'start', 'end', 'symbol', 'chrom_num', 'gene_id')
gene_info['len'] = gene_info['end'] - gene_info['start']
rownames(gene_info) = gene_info[, 'gene_id']

# harmonize the gene info ensembl ID's with the RNA-seq ensembl ID's 
info_ensembls = rownames(gene_info)
expr_ensembls = colnames(geneExpr)
shared_ensembls = set_intersection(expr_ensembls, info_ensembls)
gene_info = gene_info[rownames(gene_info) %in% shared_ensembls,]
geneExpr = geneExpr[, colnames(geneExpr) %in% shared_ensembls]

# sort by the gene ID's to get a vector and matrix that 
# are easily calculable
gene_info = gene_info[order(rownames(gene_info)),]
geneExpr = geneExpr[,order(colnames(geneExpr))]

# extract the gene lengths as a numeric vector
gene_lens = unlist(gene_info['len'])
names(gene_lens) = rownames(gene_info)

# divide each gene vector by gene length
# calculation needs to be done row-wise therefore 
# I will transpose the expression matrix before
# the division
datRates = t(geneExpr) / gene_lens

# ivide each sample vector by sample sum of rates 
# calculation needs to be done row-wise therefore 
# I will transpose the rates matrix before the 
# division (yes, another transpose). 
datTPM = t(datRates) / colSums(datRates) * 10^6

# set datExpr0 to the TPM normalized data. 
datExpr0 = log2(datTPM + 1)

# calculate good and bad samples using WGCNA method
# remove bad samples and genes
gsg = goodSamplesGenes(datExpr0, verbose = 3);

if (!gsg$allOK){
    
    if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", "))); 
    
    if (sum(!gsg$goodSamples)>0) 
        printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", "))); 
        datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes] 
}

# visualize the samples as a dendrogram to find outliers
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(16,12)
pdf(file = "results/main/sampleClustering1.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(4,6,2,2))
plot(sampleTree, main = "Sample clustering to detect outliers",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# Don't need code which is commented 
datExpr = datExpr0;
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
save(datExpr, file = "results/main/WGCNA-01-dataInput.RData")

##########################################################
######## 2. Analysis of the power parameter ##############
##########################################################

# loading the data if already pre-calculated
lnames = load(file = "results/main/WGCNA-01-dataInput.RData");
lnames

# testing a list of powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# plotting the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# plotting an horizontal line at R^2 = 0.90 
abline(h=0.90,col="red")
plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels=powers,
     cex=cex1,
     col="red")

##########################################################
######## 3. Network construction and module detection ####
##########################################################

# calculate the modules in a clockwise fashion
net = blockwiseModules(datExpr, 
                       checkMissingData = TRUE,
                       power = 16,
                       TOMType = "unsigned",
                       minModuleSize = 30,
                       maxBlockSize = 25000,
                       reassignThreshold = 0,
                       mergeCutHeight = 0.25,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "WGCNATOM",
                       verbose = 3)

# plotting the dendrogram
sizeGrWindow(15, 15)
table(net$colors)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]],
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)

# saving the results
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "results/main/WGCNA-02-networkConstruction-auto.RData")



##########################################################
############## 4. GO Enrichment Analysis #################
##########################################################
# I loaded all the required data and specifically followed section
# 3.2 of the tutorial:
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/Tutorials/anRichment-Tutorial1.pdf
# Please check the tutorial when you don't understand a step. 
# The biomaRt steps are not part of the anRichment tutorial but are 
# essential and you can look at biomaRt tutorial instead for more 
# clarification. 


# loading the data if already pre-calculated
lnames = load(file = "results/main/WGCNA-02-networkConstruction-auto.RData");
lnames


########## convert from Ensembl to Entrez ID (GO enrichment uses 
########## Entrez rather that Ensembl)
# loading the Entrez meta data using biomart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids = unlist(colnames(datExpr))
entrez_info <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", 'hgnc_symbol'),
                     mart = mart)
entrez_info = entrez_info[entrez_info['hgnc_symbol'] != "",]

# completing the conversion 
m = match(colnames(datExpr0), entrez_info$ensembl_gene_id)
go_entrez = entrez_info$entrezgene_id[m]

########## 4a. performing gene enrichment using the GO TERMS

# build the GO collection
GOcollection = buildGOcollection(organism = "human")

# perform the enrichment analysis
GOenrichment = enrichmentAnalysis(
    classLabels = moduleColors,
    identifiers = go_entrez,
    refCollection = GOcollection,
    useBackground = "given",
    threshold = 1e-4,
    thresholdType = "Bonferroni",
    getOverlapEntrez = TRUE,
    getOverlapSymbols = TRUE,
    ignoreLabels = "grey");

# display and investigate the gene enrichments
go.table.display = GOenrichment$enrichmentTable;
go.table.display$overlapGenes = shortenStrings(go.table.display$overlapGenes,
                                               maxLength = 70,
                                               split = "|");
View(go.table.display);

########## 4b. performing gene enrichment using the IMMUNE PATHWAYS
internal_cols = internalCollection(organism = 'human')

immune_enrichment = enrichmentAnalysis(classLabels = moduleColors, identifiers = go_entrez,
                                          refCollection = internal_cols,
                                          useGroups = c("ImmunePathways"),
                                          useBackground = "given",
                                          threshold = 5e-2,
                                          nBestDataSets = 3,
                                          thresholdType = "Bonferroni", 
                                          ignoreLabels = "grey");
# display and investigate the gene enrichments
imm.table.display = immune_enrichment$enrichmentTable;
imm.table.display$overlapGenes = shortenStrings(imm.table.display$overlapGenes,
                                                maxLength = 70,
                                                split = "|");
View(imm.table.display);


########## 4c. performing gene enrichment using the BLOOD PATHWAYS
blood_enrichment = enrichmentAnalysis(classLabels = moduleColors, identifiers = go_entrez,
                                       refCollection = internal_cols,
                                       useGroups = c("BloodAtlases"),
                                       useBackground = "given",
                                       threshold = 5e-2,
                                       nBestDataSets = 3,
                                       thresholdType = "Bonferroni", 
                                       ignoreLabels = "grey");

# display and investigate the gene enrichments
bl.table.display = blood_enrichment$enrichmentTable;
bl.table.display$overlapGenes = shortenStrings(bl.table.display$overlapGenes,
                                               maxLength = 70,
                                               split = "|");
View(bl.table.display);

# older way of running this method
# GOenr = GOenrichmentAnalysis(moduleColors, entrez_info$entrezgene_id[m],
# organism = "human", nBestP = 50);

########## 4d. performing gene enrichment using the MSigDB gene sets
# skipping for now but we made add it later
#msdbColl = MSigDBCollection(file = "results/refs/msigdb_v7.4.xml", organism = "human")


##########################################################
################### 5. GSEA analysis #####################
##########################################################
# GSEA requires ranks for the genes and we currently don't have 
# these types of values, for now this method will remain commented 
# until either phased in or completely removed. 
# # ## 4. Interfacing network analysis with other data such as functional
# # annotation and gene ontology
# lnames = load('results/main/WGCNA-02-networkConstruction-auto.RData')
# 
# # loading the biomart 
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# ensembl_ids = unlist(colnames(datExpr))
# entrez_info <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", 'hgnc_symbol'),
#                      mart = mart)
# entrez_info = entrez_info[entrez_info['hgnc_symbol'] != "", ]
# 
# # load the pathways file 
# gmt.file <-"results/refs/h.all.v7.4.symbols.gmt"
# pathways <- gmtPathways(gmt.file)
# 
# # convert from Ensembl to Entrez ID
# m = match(colnames(datExpr0), entrez_info$ensembl_gene_id)
# fgsea_data = datExpr0
# colnames(fgsea_data) = entrez_info$hgnc_symbol[m]
# 
# # run gsea
# fgseaRes <- fgsea(pathways = pathways, 
#                   mat = fgsea_data,
#                   minSize  = 15,
#                   maxSize  = 500)


##########################################################
################ 5. Network visualization ################
##########################################################
# ## 5. Network visualization using WGCNA functions
# dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 16);
# dissTOM[1:5,1:5]
# plotTOM = dissTOM^16;
# diag(plotTOM) = NA;
# #sizeGrWindow(12,12)
# #TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
# nSelect = 1000
# set.seed(10);
# select = sample(nGenes, size = nSelect);
# selectTOM = dissTOM[select, select];
# selectTree = hclust(as.dist(selectTOM), method = "average")
# selectColors = moduleColors[select];
# sizeGrWindow(9,9)
# plotDiss = selectTOM^16;
# diag(plotDiss) = NA;
# TOMplot(plotDiss,
#         selectTree,
#         selectColors,
#         main = "Network heatmap plot, selected genes")
# MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# MEs
# names(datTraits)
# MET = orderMEs(MEs)
# sizeGrWindow(5,7.5);
# par(cex = 0.9)
# plotEigengeneNetworks(MET, "", 
#                       marDendro = c(0,4,1,2),
#                       marHeatmap = c(3,4,1,2),
#                       cex.lab = 0.8,
#                       xLabelsAngle = 90)