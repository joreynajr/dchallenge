rule process_split_ols_ids:
    input:
        'results/notebooks/gwas_catalog_disease_traits/splits/split.ols_id.{split_no}.txt'
    output:
        'results/notebooks/gwas_catalog_disease_traits/processed/split.ols_id.{split_no}.txt'
    log: 
        'results/notebooks/gwas_catalog_disease_traits/logs/process_split_ols_ids.{split_no}.log'
    resources:
        mem_mb = 6000,
        nodes = 1,
        ppn = 1,
    shell:
        r'''
            {config[hichip_db_python]} \
                        workflow/notebooks/exploratory/GWAS_Catalog_Disease_Traits.Search_API.py \
                                    {input} {output} >> {log} 2>&1
        '''
        

# liftover GWAS summary statistics to GRCh38 when necessary.
rule liftover_grch37_to_38_gwas_summary_stats: 
    input:
        ss = 'results/main/coloc/Data/T1D_GWAS/{gwas_source}/GRCh37/GWAS_input_colocalization_pval_lt_5eMinus8.txt',
        chain = 'results/refs/ucsc/hg19ToHg38.over.chain.gz'
    output:
        'results/main/coloc/Data/T1D_GWAS/{gwas_source}/GRCh38/GWAS_input_colocalization_pval_lt_5eMinus8.txt',
    log:
        'results/main/coloc/Data/T1D_GWAS/{gwas_source}/logs/rule.liftover_grch37_gwas_summary_stats.{gwas_source}.log'
    wildcard_constraints:
        gwas_source = "\w"
    shell:
        r'''
            # convert the SNP file into bed format 
            echo # convert the SNP file into bed format >> {log}
            lift_input="{input.ss}.input.txt"
            cat {input.ss} | \
                    sed '1d' | \
                    awk 'BEGIN{{OFS="\t"}}; {{print $1, ($2 - 1), $2, $3, $4, $5}}' > \
                    $lift_input

            # liftover from hg37 to hg38
            echo # liftover from hg37 to hg38 >> {log}
            lifted="{input.ss}.mapped.bed"
            unmapped="{input.ss}.unmapped.bed"
            /home/jreyna/software/UCSC_Browser_Tools/liftOver \
                                -bedPlus=3 \
                                -tab \
                                $lift_input \
                                {input.chain} \
                                $lifted \
                                $unmapped

            # make the final input file
            echo # make the final input file >> {log}
            # add the header back and convert the 
            head -n 1 {input.ss} | awk 'BEGIN{{OFS="\t"}}; {{print $0}}' > {output}
            awk 'BEGIN{{OFS=" "}}; {{print $1, $3, $4, $5, $6}}' $lifted >> {output}

            # remove intermediate files
            rm $lift_input $lifted $unmapped
        '''
