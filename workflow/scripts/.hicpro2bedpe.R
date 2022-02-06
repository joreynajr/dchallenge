library(HiCcompare)
args = commandArgs(trailingOnly=TRUE)

# read in data
mat <- read.table(args[1])
bed <- read.table(args[2])

# convert to BEDPE
dat <- hicpro2bedpe(mat, bed)

write.table(dat, file=args[3], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
