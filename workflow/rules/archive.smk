#########################################################################################
# Rules to process eQTL Catalogue with eQTL FDR filtering
# Temporarily Deprecated: Requires synchronization/update with previous rules
#########################################################################################
# Processing the eQTL files to extract only significant eQTLs. 
# Currently this is not needed an a more streamline approach is available
# by using the already formatted colocalization input from rule. 
# processing_eqtl_catalog_with_complete_fields. This rule will
# be deprected for the time being but is required for rule
# calculate_number_of_eqtls_post_filtering:
#rule processing_eqtl_catalog_with_complete_fields_and_FDR_filtering: # (status: running)
#    input:
#        rules.add_missing_cols.output
#    output:
#        'results/main/eqtl/{eqtl_source}/ge/{eqtl_source}_ge_{ge_source}.all.complete_fields.input.tsv.gz'
#    params:
#        interm='results/main/eqtl/{eqtl_source}/ge/{eqtl_source}_ge_{ge_source}.all.complete_fields.input.tsv'
#    log: 
#        'results/main/coloc/Results/eQTL_Catalogue_FDR/logs/parse_eqtl_catalog_with_complete_fields.{eqtl_source}.{ge_source}.log'
#    resources:
#        mem_mb = 12000,
#        nodes = 1,
#        ppn = 1,
#    shell:
#        r'''
#            echo "zcat + sed + awk processing" >> {log}
#            zcat {input} | \
#                 sed '1d' | \
#                 awk 'BEGIN{{OFS="	"}}; {{if($21 < 0.05){{print "chr"$2, $3, $17, $20, $10, $9, $21}}}}' > \
#                 {params.interm} 2>> {log}
#            {config[bgzip]} {params.interm} 2>> {log}
#        '''
#
## Count the number of eqtl after filtering with FDR
#rule calculate_number_of_eqtls_post_filtering: # (Status: running)
#    input:
#        eqtl = rules.processing_eqtl_catalog_with_complete_fields.output
#    output:
#        'results/main/eqtl/{eqtl_source}/ge/{eqtl_source}_ge_{ge_source}.all.postfilter.num_eqtls.txt'
#    log: 
#        'results/main/eqtl/logs/calculate_number_of_eqtls_post_filtering.{eqtl_source}.{ge_source}.log'
#    resources:
#        mem_mb = 4000,
#        nodes = 1,
#        ppn = 1,
#    shell:
#        r'''
#            zcat {input} | wc -l  > {output} 2> {log}
#        '''
#
## run colocalization with eqtls which are significant (filter on FDR < 0.05)
#rule run_colocalization_eqtl_catalog_with_eqtl_fdr_filtering:
#    input:
#        gwas = 'results/main/coloc/Data/T1D_GWAS/{gwas_source}/GRCh38/GWAS_input_colocalization_pval_lt_5eMinus8.txt',
#        eqtl = rules.processing_eqtl_catalog_with_complete_fields.output
#    output:
#        outdir = directory('results/main/coloc/Results/eQTL_Catalogue_FDR/{gwas_source}/{eqtl_source}/{ge_source}/')
#    log: 
#        'results/main/coloc/Results/eQTL_Catalogue_FDR/logs/run_colocalization_eqtl_catalog.{gwas_source}.{eqtl_source}.{ge_source}.log'
#    params:
#        chr = 1,
#        pos = 2,
#        gene_name = 3,
#        dist = 4,
#        slope = 5,
#        pvalue = 6,
#        fdr = 7,
#        header = 'FALSE'
#    resources:
#        mem_mb = 48000,
#        nodes = 1,
#        ppn = 1,
#    shell:
#        r"""
#            # make the log directory
#            mkdir -p $(dirname {log})
#
#            # set the tempdir 
#            TMPDIR=/scratch
#
#            # add the corresponding R for colocalization analyses
#            PATH=/mnt/BioApps/R/3.6.1/bin/:$PATH
#
#            # run colocalization command
#            # omitting --eqtl-FDR {params.fdr} since eQTL Catalog doesn't have FDR/qvalue
#            # Omitting --eqtl-header {params.header} since header was removed previously
#            Rscript scripts/coloc/Colocalization_Analysis_GWAS_Script_Generalized.R \
#                                --eqtl-chr {params.chr} \
#                                --eqtl-pos {params.pos} \
#                                --eqtl-geneName {params.gene_name} \
#                                --eqtl-dist {params.dist} \
#                                --eqtl-slope {params.slope} \
#                                --eqtl-pvalue {params.pvalue} \
#                                --eqtl-FDR {params.fdr} \
#                                {input.gwas} \
#                                {input.eqtl} \
#                                {output.outdir} >> {log} 2>&1
#        """


##########################################################################################
## Temporarily Deprecated: Underdeveloped ImmuNexUT Rules
##########################################################################################
#
## This file doesn not contain a chr and position field, rather it has a 
## sid filed form which you can generated these two.
#rule parse_immunexut_eqtl: #(Status: developing)
#    input:
#        'results/main/eqtl/ImmuNexUT/raw/E-GEAD-420/{ge_source}_nominal.txt'
#    output:
#        'results/main/eqtl/ImmuNexUT/ge/ImmuNexUT_ge_{ge_source}.all.immunenexut_pre.tsv'
#    log: 
#        'results/main/coloc/Results/ImmuNexUT_FDR/logs/parse_immunexut_eqtl.ImmuNexUT.{ge_source}.log'
#    shell:
#        r'''
#            echo "chr	pos	genename	dist	slope	pval" > {output}
#            cat {input} | sed '1d' | \
#                awk 'BEGIN{{OFS="	"}} function abs(v) {{return v < 0 ? -v : v}} \
#                    {{dist=abs($4 - $10); print $9, $10, $1, dist, $13, $12}}' >> \
#                {output} 2> {log}
#        '''
#
#
## This file doesn not contain a chr and position field, rather it has a 
## sid field form which you can generated these two.
#rule calculate_fdr_col_immunexut: # (Status: developing)
#    input:
#        eqtl = rules.parse_immunexut_eqtl.output
#    output:
#        'results/main/eqtl/{eqtl_source}/ge/{eqtl_source}_ge_{ge_source}.all.immunexut.fdr.txt'
#    log: 
#        'results/main/eqtl/logs/calculate_fdr_col_immunexut.{eqtl_source}.{ge_source}.log'
#    params:
#        pval = 6
#    resources:
#        mem_mb = 24000,
#        nodes = 1,
#        ppn = 1,
#    shell:
#        r'''
#            /mnt/BioHome/jreyna/software/anaconda3/envs/hic_tls/bin/python 
#                        workflow/scripts/eqtl/Calc_FDR_Col.py \
#                        --in-header \
#                        --out-header \
#                        {input} \
#                        {output} \
#                        {params.pval} >> {log} 2>&1
#        '''
#
#
## After running calculate_dist_col this rule will add the dist column to the tsv.gz data
#rule add_missing_cols_immunexut: #(Status: developing)
#    input:
#        eqtl = rules.parse_immunexut_eqtl.output,
#        fdr = rules.calculate_fdr_col_immunexut.output
#    params:
#        interm = 'results/main/eqtl/{eqtl_source}/ge/{eqtl_source}_ge_{ge_source}.all.immunexut.complete_fields.tsv'
#    output:
#        'results/main/eqtl/{eqtl_source}/ge/{eqtl_source}_ge_{ge_source}.all.immunexut.complete_fields.tsv.gz'
#    log: 
#        'results/main/eqtl/logs/add_missing_cols_immunexut.{eqtl_source}.{ge_source}.log'
#    resources:
#        mem_mb = 32000,
#        nodes = 1,
#        ppn = 1,
#    shell:
#        r'''
#            echo "running paste" >> {log} 2>&1
#            paste {input.eqtl} {input.fdr} | tr -d '\r' > {params.interm} 2>> {log}
#
#            echo "running bgzip" >> {log} 2>&1
#            {config[bgzip]} {params.interm} 2>> {log}
#        '''
#
#
#
## This file doesn not contain a chr and position field, rather it has a 
## sid filed form which you can generated these two.
#rule processing_immunexut_with_complete_fields: # (status: developing)
#    input:
#        rules.add_missing_cols_immunexut.output
#    output:
#        'results/main/eqtl/{eqtl_source}/ge/{eqtl_source}_ge_{ge_source}.all.immunexut.complete_fields.input.tsv.gz'
#    params:
#        interm='results/main/eqtl/{eqtl_source}/ge/{eqtl_source}_ge_{ge_source}.all.immunexut.complete_fields.input.tsv'
#    log: 
#        'results/main/coloc/Results/ImmuNexUT/logs/parse_eqtl_catalog_with_complete_fields.{eqtl_source}.{ge_source}.log'
#    resources:
#        mem_mb = 12000,
#        nodes = 1,
#        ppn = 1,
#    shell:
#        r'''
#            #awk 'BEGIN{{OFS="	"}}; {{if($21 < 0.05){{print "chr"$1, $2, $3, $4, $5, $6, $$7}}}}' > \
#            echo "zcat + sed + awk processing" >> {log}
#            zcat {input} | \
#                 sed '1d' | \
#                 awk 'BEGIN{{OFS="	"}}; {{print $1, $2, $3, $4, $5, $6, $$7}}' > \
#                 {params.interm} 2>> {log}
#            {config[bgzip]} {params.interm} 2>> {log}
#        '''
#
#
#use rule run_colocalization_eqtl_catalog as run_colocalization_eqtl_catalog_immunexut with: # (Status: developing)
#    input:
#        gwas = 'results/main/coloc/Data/T1D_GWAS/{gwas_source}/GRCh38/GWAS_input_colocalization_pval_lt_5eMinus8.txt',
#        eqtl = rules.processing_immunexut_with_complete_fields.output
#    output:
#        outdir = directory('results/main/coloc/Results/ImmuNexUT/{gwas_source}/{eqtl_source}/{ge_source}/')
#    log: 
#        'results/main/coloc/Results/ImmuNexUT/logs/run_colocalization_eqtl_catalog.{gwas_source}.{eqtl_source}.{ge_source}.log'
#    params:
#        chr = 1,
#        pos = 2,
#        gene_name = 3,
#        dist = 4,
#        slope = 5,
#        pvalue = 6,
#        fdr = 'NULL',
#        header = 'FALSE'
#
#
#use rule run_colocalization_eqtl_catalog_with_fdr_filtering as run_colocalization_eqtl_catalog_with_fdr_filtering_immunexut with:
#    input:
#        gwas = 'results/main/coloc/Data/T1D_GWAS/{gwas_source}/GRCh38/GWAS_input_colocalization_pval_lt_5eMinus8.txt',
#        eqtl = rules.processing_immunexut_with_complete_fields.output
#    output:
#        outdir = directory('results/main/coloc/Results/ImmuNexUT_FDR/{gwas_source}/{eqtl_source}/{ge_source}/')
#    log: 
#        'results/main/coloc/Results/ImmuNexUT_FDR/logs/run_colocalization_eqtl_catalog.{gwas_source}.{eqtl_source}.{ge_source}.log'
#
#
## This file doesn not contain a chr and position field, rather it has a 
## sid filed form which you can generated these two.
#rule parse_immunexut_eqtl: #(Status: developing)
#    input:
#        'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.txt.gz'
#    output:
#        'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.input.txt'
#    params:
#        unzipped = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.txt',
#        interm1 = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.interm1.txt',
#        interm2 = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.interm2.txt'
#    log: 
#        'results/main/coloc/Results/ImmuNexUT/logs/parse_mu_et_al_eqtl.{eqtl_source}.{ge_source}.log'
#    shell:
#        r'''
#            zcat {input} > {params.unzipped}
#             cut -f 6 {params.unzipped} | \
#                 sed '1d' | \
#                 awk 'BEGIN{{FS=":"; OFS=" "; print "chr   pos"}} {{print "chr"$1,$2}}' > {params.interm1}
#            paste {params.unzipped} {params.interm1} > {params.interm2}
#            awk 'BEGIN{{OFS="	"}} {{print $13, $14, $1, $7, $9, $11, $12}}' {params.interm2} > {output}
#            rm {params.unzipped} {params.interm1} {params.interm2}
#        '''

#########################################################################################
# Semi-permanent Deprecated: Using eQTL Catalogue and ImmuNexUT instead
#########################################################################################
## This file doesn not contain a chr and position field, rather it has a 
## sid filed form which you can generated these two.
#rule parse_mu_et_al_eqtl: #(Status: developing)
#    input:
#        'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.txt.gz'
#    output:
#        'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.input.txt'
#    params:
#        unzipped = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.txt',
#        interm1 = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.interm1.txt',
#        interm2 = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.interm2.txt'
#    log: 
#        'results/main/coloc/Results/Mu_et_al/logs/parse_mu_et_al_eqtl.{eqtl_source}.{ge_source}.log'
#    shell:
#        r'''
#            zcat {input} > {params.unzipped}
#             cut -f 6 {params.unzipped} | \
#                 sed '1d' | \
#                 awk 'BEGIN{{FS=":"; OFS=" "; print "chr   pos"}} {{print "chr"$1,$2}}' > {params.interm1}
#            paste {params.unzipped} {params.interm1} > {params.interm2}
#            awk 'BEGIN{{OFS="	"}} {{print $13, $14, $1, $7, $9, $11, $12}}' {params.interm2} > {output}
#            rm {params.unzipped} {params.interm1} {params.interm2}
#        '''
#
#
##rules.parse_mu_et_al_eqtl.output:
#rule run_colocalization_mu_et_al: # (Status: developing)
#    input:
#        gwas = 'results/main/coloc/Data/T1D_GWAS/{gwas_source}/GRCh38/GWAS_input_colocalization_pval_lt_5eMinus8.txt',
#        eqtl = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.input.txt'
#    output:
#        outdir = directory('results/main/coloc/Results/Mu_et_al/{gwas_source}/{eqtl_source}/{ge_source}/')
#    log: 
#        'results/main/coloc/Results/Mu_et_al/logs/run_colocalization_mu_et_al_example.{gwas_source}.{eqtl_source}.{ge_source}.log'
#    params:
#        chr = 1,
#        pos = 2,
#        gene_name = 3,
#        dist = 4,
#        slope = 5,
#        pvalue = 6,
#        fdr = 7,
#        header = 'TRUE'
#    resources:
#        mem_mb = 20000,
#        nodes = 1,
#        ppn = 1,
#        #PBS -j eo
#        #PBS -V
#    shell:
#        r"""
#            # make the log directory
#            mkdir -p $(dirname {log})
#
#            # set the tempdir 
#            TMPDIR=/scratch
#
#            # add the corresponding R for colocalization analyses
#            PATH=/mnt/BioApps/R/3.6.1/bin/:$PATH
#
#            # run colocalization command
#            Rscript scripts/coloc/Colocalization_Analysis_GWAS_Script_Generalized.R \
#                                --eqtl-chr {params.chr} \
#                                --eqtl-pos {params.pos} \
#                                --eqtl-geneName {params.gene_name} \
#                                --eqtl-dist {params.dist} \
#                                --eqtl-slope {params.slope} \
#                                --eqtl-pvalue {params.pvalue} \
#                                --eqtl-FDR {params.fdr} \
#                                --eqtl-header {params.header} \
#                                {input.gwas} \
#                                {input.eqtl} \
#                                {output.outdir} >> {log} 2>&1
#        """
