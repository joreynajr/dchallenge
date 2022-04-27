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
