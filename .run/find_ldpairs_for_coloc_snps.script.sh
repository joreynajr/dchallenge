fns=""
ref_genome="GRCh37"

echo "Is it running?"
tracker=".run/coloc_ldpairs.sh.tracker.txt"
rscript="/share/apps/R/3.4.3/lib64/R/bin/Rscript"

while IFN= read -r line;
do

    #fn = 'results/main/coloc/eQTL_Catalogue/{gwas_source}/{eqtl_source}/{ge_source}/ldpairs/coloc_ld_snps.txt'
    #fn="results/main/coloc/eQTL_Catalogue/${run_info[0]}/${ref_genome}/${run_info[1]}/${run_info[2]}/ldpairs/coloc_ld_snps.txt"
    #fns="$fn "
    #echo $fn

    # extracting values 
    # run_info = GWAS, eQTL DB, eQTL origin
    run_info=($line)
    gwas_source=${run_info[0]}
    eqtl_source=${run_info[1]}
    ge_source=${run_info[2]}

    msg="Processing ${gwas_source} ${eqtl_source} ${ge_source}"
    echo "*******************************************************************************"
    echo $msg
    echo "*******************************************************************************"

    # setting input
    snp_file="results/main/coloc/Results/eQTL_Catalogue/${gwas_source}/${eqtl_source}/${ge_source}/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed"
    onekg_dir="/mnt/BioAdHoc/Groups/vd-vijay/Ariel/R24_new/LD-snps-1000G/DupsRemoved/"
    population_dir="/mnt/BioAdHoc/Groups/vd-vijay/Ariel/R24_new/LD-snps/lists-pops/"
    snpinfo_dir="/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2020_IQTL_HiChIP/Data/SNPInfo/SNPInfo_merged_tables/"

    # setting params
    workdir="results/main/coloc/eQTL_Catalogue/${gwas_source}/${ref_genome}/${eqtl_source}/${ge_source}/ldpairs/"
    chr_col=1
    pos_col=2

    # setting output
    ld="results/main/coloc/eQTL_Catalogue/${gwas_source}/${ref_genome}/${eqtl_source}/${ge_source}/ldpairs/coloc_ld_snps.txt"
    log="results/main/coloc/eQTL_Catalogue/${gwas_source}/${ref_genome}/${eqtl_source}/${ge_source}/ldpairs/coloc_ld_snps.log"
    outfn="results/main/coloc/eQTL_Catalogue/${gwas_source}/${ref_genome}/${eqtl_source}/${ge_source}/ldpairs/Out_Merge_LD.txt"

    wildcards="${gwas_source}/${ref_genome}/${eqtl_source}/${ge_source}/"
    if [[ ! -f "$snp_file" ]];
    then
        echo "SKIPPING: $snp_file is not available now, skipping $wildcards." | tee -a $tracker
        continue
    elif [[ -f "outfn" ]];
    then
        echo "SKIPPING: $outfn was already generated, skipping $wildcards." | tee -a $tracker
        continue
    else
        echo "RUNNING: $ld, running $wildcards." | tee -a $tracker
    fi

    # prining the command for logging
    cmd_str="$rscript workflow/scripts/ldpairs/ldpair_with_plink.R \n\t
                --snp-file ${snp_file} \n\t
                --onekg-dir ${onekg_dir} \n\t
                --population-dir ${population_dir} \n\t
                --snpinfo-dir ${snpinfo_dir} \n\t
                --header \n\t
                --chr-col ${chr_col} \n\t
                --chr-prefix \n\t
                --pos-col ${pos_col} \n\t
                --workdir ${workdir}"

    #echo "Running the command:"
    #echo -e $cmd_str

    # running the command
    $rscript workflow/scripts/ldpairs/ldpair_with_plink.R \
                --snp-file ${snp_file} \
                --onekg-dir ${onekg_dir} \
                --population-dir ${population_dir} \
                --snpinfo-dir ${snpinfo_dir} \
                --header \
                --chr-col ${chr_col} \
                --chr-prefix \
                --pos-col ${pos_col} \
                --workdir ${workdir}

    old_fn="${workdir}/Out_Merge_LD.txt"
    mv $old_fn ${ld}

done < <(cat config/coloc_samplesheets/t1d.samplesheet.tsv | sed '1d' | grep -v "^#")

#echo $fns
#snakemake --profile workflow/profiles/local $@ $fns
#snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
