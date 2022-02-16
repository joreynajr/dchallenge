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
        
