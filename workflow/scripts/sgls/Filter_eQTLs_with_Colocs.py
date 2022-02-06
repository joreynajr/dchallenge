#!/usr/bin/env python
# coding: utf-8
import sys
import gzip

# test
#eqtl = 'results/main/eqtl/Schmiedel_2018/ge/Schmiedel_2018_ge_monocyte_naive.all.complete_fields.tsv.gz'
#coloc = 'results/main/coloc/Results/Colocalization_SMKN/T1D_34012112_Gaulton/Schmiedel_2018/monocyte_naive/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed'
#outfn = 'results/main/sgls/T1D_34012112_Gaulton/Schmiedel_2018/monocyte_naive/monocyte_naive/eqtls.coloc_filtered.tsv.gz'

eqtl = sys.argv[1]
coloc = sys.argv[2]
outfn = sys.argv[3]

# find all the coloc ID's composed of the chromosome:position
coloc_ids = []
with open(coloc) as fr:
    
    # skip the first line
    next(fr)
    
    # process the remaining lines
    for line in fr:
        info = line.split('\t')
        
        chrom = info[0].replace('chr', '')
        pos = info[1]
        idd = '{}:{}'.format(chrom, pos)
        coloc_ids.append(idd)
coloc_ids = set(coloc_ids)

# filter/extract eQTLs SNPs within the list of coloc_ids or fdr < 0.05
with gzip.open(eqtl, mode='rt') as fr, gzip.open(outfn, mode='wb') as fw:

    # get and write the header
    header = next(fr)
    fw.write(header.encode())

    # parse the remaining lines
    for line in fr:
        info = line.split()
        
        idd = '{}:{}'.format(*info[1:3])
        fdr = float(info[-1])
        
        if idd in coloc_ids or fdr < 0.05:
            fw.write(line.encode())
        
