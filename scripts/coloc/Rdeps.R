# Install the R packages required for the colocalization analysis

# libraries installed with CRAN
install.packages("data.table")
install.packages("stringr")
install.packages("coloc")
install.packages("ggplot2")
install.packages("ggplotify")

# libraries installed with Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
library("BiocManager")
BiocManager::install("snpStats")
