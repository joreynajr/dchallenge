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


#use rule finemapped_ldpairs as coloc_ldpairs with:
rule coloc_ldpairs:
    input:
        snp_file = 'results/main/coloc/Results/eQTL_Catalogue/{gwas_source}/{eqtl_source}/{ge_source}/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed',
        onekg_dir = '/mnt/BioAdHoc/Groups/vd-vijay/Ariel/R24_new/LD-snps-1000G/DupsRemoved/',
        population_dir = '/mnt/BioAdHoc/Groups/vd-vijay/Ariel/R24_new/LD-snps/lists-pops/',
        snpinfo_dir = '/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2020_IQTL_HiChIP/Data/SNPInfo/SNPInfo_merged_tables/'
    params:
        workdir = 'results/main/coloc/eQTL_Catalogue/{gwas_source}/{ref_genome}/{eqtl_source}/{ge_source}/ldpairs/',
        chr_col = 1,
        pos_col = 2
    resources:
        mem_mb = 20000
    output:
        ld = 'results/main/coloc/eQTL_Catalogue/{gwas_source}/{ref_genome}/{eqtl_source}/{ge_source}/ldpairs/coloc_ld_snps.txt'
    log:
        'results/main/coloc/eQTL_Catalogue/{gwas_source}/{ref_genome}/{eqtl_source}/{ge_source}/ldpairs/coloc_ld_snps.log'
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
