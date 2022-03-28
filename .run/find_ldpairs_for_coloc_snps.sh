fns=""
ref_genome="GRCh37"

while IFN= read -r line;
do

    # extracting values 
    # run_info = GWAS, eQTL DB, eQTL origin
    run_info=($line)
    gwas_source=${run_info[0]}
    eqtl_source=${run_info[1]}
    ge_source=${run_info[2]}

    # getting the current LD target
    ld="results/main/GRCh37/coloc/eQTL_Catalogue/${gwas_source}/${eqtl_source}/${ge_source}/ldpairs/coloc_ld_snps.txt"

    # calculate LD if the colocalization data has previously been processed
    coloc="results/main/coloc/Results/eQTL_Catalogue/${gwas_source}/${eqtl_source}/${ge_source}/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed"
    if [ -f "$coloc" ];
    then
        fns+="$ld "
        #break
    fi

done < <(cat config/coloc_samplesheets/t1d.samplesheet.tsv | sed '1d' | grep -v "^#")

echo $fns
#snakemake --profile workflow/profiles/local $@ $fns
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
