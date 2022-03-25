rule run_finemapping:
    input:
        'results/main/finemapping/{gwas_source}/{ref_genome}/raw/GWAS_RawData.finemap.input.tsv'
    output:
        'results/main/finemapping/{gwas_source}/{ref_genome}/offset_{offset}/Summary/sss/FINAL_top_snp_credible_set.txt'
    params:
        chrcol = 1,
        poscol = 2,
        pvalcol = 3,
        afcol = 4,
        betacol = 5,
        secol = 6 
    resources:
        mem_mb = 10000,
        nodes = 1,
        ppn = 8
    log:
        'results/main/finemapping/logs/rule_run_finemapping.{gwas_source}.{ref_genome}.offset_{offset}.log'
    shell:
        r'''
            # make a config file
            mkdir -p {output}
            abs_config=$(readlink -f "{output}/config.txt")
            abs_input=$(readlink -f {input})
            abs_output=$(readlink -f {output})
            cat workflow/scripts/finemapping/Finemap_SMK/configfile.template.txt | \
                        sed "s+{{gwas}}+$abs_input+" | \
                        sed "s+{{outdir}}+$abs_output+" | \
                        sed "s+{{offset}}+{wildcards.offset}+" | \
                        sed "s+{{chrcol}}+{params.chrcol}+" | \
                        sed "s+{{poscol}}+{params.poscol}+" | \
                        sed "s+{{pvalcol}}+{params.pvalcol}+" | \
                        sed "s+{{afcol}}+{params.afcol}+" | \
                        sed "s+{{betacol}}+{params.betacol}+" | \
                        sed "s+{{secol}}+{params.secol}+"  > $abs_config

            # run the finemapping
            cd workflow/scripts/finemapping/Finemap_GWAS_Joaquin/
            bash Finemap_Script.sh -C $abs_config
        '''

#rule dev_ldpairs:
#    input:
#        snp_file = 'results/main/finemapping/T1D_32005708/GRCh37/offset_1000000/Summary/sss/FINAL_top_snp_credible_set.txt',
#        onekg_dir = '/mnt/BioAdHoc/Groups/vd-vijay/Ariel/R24_new/LD-snps-1000G/DupsRemoved/',
#        population_dir = '/mnt/BioAdHoc/Groups/vd-vijay/Ariel/R24_new/LD-snps/lists-pops/',
#        snpinfo_dir = '/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2020_IQTL_HiChIP/Data/SNPInfo/SNPInfo_merged_tables/'
#    params:
#        workdir ='tests/test_ldpairing/',
#        chr_col = 5,
#        pos_col = 6
#    output:
#        ld = 'tests/test_ldpairing/ld_snps.txt'
#    log:
#        'tests/test_ldpairing/test.log'
#    shell:
#        r"""
#            Rscript workflow/scripts/ldpairs/ldpair_with_plink.R \
#                        --snp-file {input.snp_file} \
#                        --onekg-dir {input.onekg_dir}/ \
#                        --population-dir {input.population_dir}/ \
#                        --snpinfo-dir {input.snpinfo_dir}/ \
#                        --header \
#                        --chr-col {params.chr_col} \
#                        --pos-col {params.pos_col} \
#                        --workdir {params.workdir}/
#            old_fn="{params.workdir}/Out_Merge_LD.txt"
#            mv $old_fn {output}
#        """

# Performing LD Pair analysis for finemapped data
rule finemapped_ldpairs:
    input:
        snp_file = 'results/main/finemapping/{gwas_source}/{ref_genome}/offset_{offset}/Summary/sss/FINAL_top_snp_credible_set.txt',
        onekg_dir = '/mnt/BioAdHoc/Groups/vd-vijay/Ariel/R24_new/LD-snps-1000G/DupsRemoved/',
        population_dir = '/mnt/BioAdHoc/Groups/vd-vijay/Ariel/R24_new/LD-snps/lists-pops/',
        snpinfo_dir = '/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2020_IQTL_HiChIP/Data/SNPInfo/SNPInfo_merged_tables/'
    params:
        workdir ='results/main/finemapping/{gwas_source}/{ref_genome}/offset_{offset}/ldpairs/',
        chr_col = 5,
        pos_col = 6
    resources:
        mem_mb = 20000
    output:
        ld = 'results/main/ldpairs/finemapping/{gwas_source}/{ref_genome}/offset_{offset}/finemapping_ld_snps.txt'
    wildcard_constraints:
        offset = '[0-9]*'
    log:
        'results/main/finemapping/{gwas_source}/{ref_genome}/offset_{offset}/ldpairs/finemapping_ld_snps.log'
    shell:
        r"""
            Rscript workflow/scripts/ldpairs/ldpair_with_plink.R \
                        --snp-file {input.snp_file} \
                        --onekg-dir {input.onekg_dir}/ \
                        --population-dir {input.population_dir}/ \
                        --snpinfo-dir {input.snpinfo_dir}/ \
                        --header \
                        --chr-col {params.chr_col} \
                        --pos-col {params.pos_col} \
                        --workdir {params.workdir}/
            old_fn="{params.workdir}/Out_Merge_LD.txt"
            mv $old_fn {output}
        """
