fns=""
for fn in $(ls results/main/coloc/Results/eQTL_Catalogue/*/*/*/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed);
do
    # get the new file name
    stem="results/main/GRCh37/coloc/Results/eQTL_Catalogue/"
    rest=$(echo $fn | cut -d "/" -f 6-)
    new_fn="$stem/$rest"

    #echo $fn
    #echo $new_fn
    fns+="$new_fn "
done

#snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
snakemake --profile workflow/profiles/local/ $@ $fns
