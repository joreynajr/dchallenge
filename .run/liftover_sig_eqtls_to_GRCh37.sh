fns=""
eqtl_fns=$(ls results/main/sgls/*/*/*/eqtls.coloc_filtered.tsv.gz | grep -v sgls_20220222)
for eqtl_fn in $eqtl_fns;
do

    # extract the run information
    run_info=($(echo $eqtl_fn | cut -d "/" -f 4-6 | sed "s/\// /g"))

    # set up the colocalization path
    coloc_fn="results/main/coloc/Results/eQTL_Catalogue/${run_info[0]}/${run_info[1]}/${run_info[2]}/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed"

    # run when you have both the eQTL and colocalization data
    if [ -f "$eqtl_fn" ] && [ -f "$coloc_fn" ];
    then
        stem="results/main/GRCh37/sgls/"
        rest=$(echo $eqtl_fn | cut -d "/" -f 4-)
        new_fn="$stem/$rest"
        fns+="$new_fn "
    fi
done

#echo $fns
#snakemake --profile workflow/profiles/local/ $@ $fns
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
