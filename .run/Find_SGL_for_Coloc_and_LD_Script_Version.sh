fns=""
while IFS= read -r line;
do
    # get the params
    run_info=($line)
    gwas_source="${run_info[0]}"
    eqtl_source="${run_info[1]}"
    ge_source="${run_info[2]}"
    loop_source="${run_info[3]}"

    #echo "creating the following: ${run_info[0]}/GRCh37/${run_info[1]}/${run_info[2]}/${run_info[3]}/"

    # get the dependency files
    eqtl="results/main/GRCh37/sgls/${gwas_source}/${eqtl_source}/${ge_source}/eqtls.coloc_filtered.tsv.gz"
    coloc="results/main/GRCh37/coloc/eQTL_Catalogue/${gwas_source}/${eqtl_source}/${ge_source}/ldpairs/coloc_ld_snps.txt"
    loops="results/main/h3k27ac_hichip/${loop_source}/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed"
    output_fn="results/main/GRCh37/sgls/ldpairs/${gwas_source}/${eqtl_source}/${ge_source}/${loop_source}/script_version/"

    # add to the runnable files if all deps are met
    if [ -f "$eqtl" ] && [ -f "$coloc" ] && [ -f "$loops" ];
    then
        echo "PROCESSING: $output_fn"
        fns+="$output_fn "
    else
        echo "WILL NOT PROCESS: $output_fn"
        echo "WARNING: Missing dependencies $eqtl $coloc $loops"
    fi

done < <(cat config/sgl_samplesheets/sgl.samplesheet.t1d_working.tsv | sed '1d' | grep -v "^#" | head -n 1000)

echo "fns: $fns"
snakemake --profile workflow/profiles/local/ $@ $fns
#snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
