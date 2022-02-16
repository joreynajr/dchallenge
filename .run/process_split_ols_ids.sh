outfns=""
for i in $(seq 0 36);
do
    outfn="results/notebooks/gwas_catalog_disease_traits/processed/split.ols_id.$i.txt"
    outfns+="$outfn "
done

echo $outfns
snakemake --profile workflow/profiles/pbs-torque/ $@ $outfns
