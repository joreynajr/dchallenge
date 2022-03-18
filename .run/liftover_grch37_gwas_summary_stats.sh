fns=""
for gwas in $(ls -d results/main/coloc/Data/T1D_GWAS/* | head -n 1);
do
    gwas_source=$(basename $gwas)
    fn="results/main/coloc/Data/T1D_GWAS/${gwas_source}/GRCh38/GWAS_input_colocalization_pval_lt_5eMinus8.txt"
    fns+="$fn "
done

snakemake --profile workflow/profiles/local $@ $fns
