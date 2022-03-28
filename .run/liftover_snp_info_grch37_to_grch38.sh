outfns=""
chrs=($(seq 1 22) "X")
for chr_num in ${chrs[@]};
do
    outfns+="results/refs/coloc_snps_meta/grch38/SNPInfo/snpinfo_chr${chr_num}.txt "
done
snakemake --profile workflow/profiles/pbs-torque/ $@ $outfns
