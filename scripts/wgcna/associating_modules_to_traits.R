library("WGCNA")
lnames1 = load("results/main/WGCNA-01-dataInput.RData") 
lnames2 = load("results/main/WGCNA-02-networkConstruction-auto.RData")

# Load the traits data
fn = "results/refs/datasets/GSE123658/GSE123658_sample_trait_data.tsv"
datTraits = read.table(fn, header = T, sep='\t')
datTraits = datTraits[c(1,4,7,8)]

# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
rownames(MEs) = gsub("X", "", rownames(MEs))
rownames(datExpr) = gsub("X", "", rownames(datExpr))
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# Reorder the rows so MEs and datTraits
reorder_idxs = match(rownames(MEs), datTraits$subject_id)
datTraits = datTraits[reorder_idxs,2:4]

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
