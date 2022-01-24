## This file doesn not contain a chr and position field, rather it has a 
## sid filed form which you can generated these two.
#rule parse_mu_et_al_eqtl_example:
#    input:
#        'results/main/coloc/Data/eqtl_sqtl_summ_stats/BLUEPRINT_eQTL/Monocyte.txt.gz'
#    output:
#        'results/main/coloc/Data/eqtl_sqtl_summ_stats/BLUEPRINT_eQTL/Monocyte.input.txt'
#    params:
#        unzipped = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/BLUEPRINT_eQTL/Monocyte.txt',
#        interm1 = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/BLUEPRINT_eQTL/Monocyte.interm1.txt',
#        interm2 = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/BLUEPRINT_eQTL/Monocyte.interm2.txt'
#    shell:
#        r'''
#            zcat {input} > {params.unzipped}
#             cut -f 6 {params.unzipped} | \
#                 sed '1d' | \
#                 awk 'BEGIN{{FS=":"; OFS=" "; print "chr   pos"}} {{print $1,$2}}' > {params.interm1}
#            paste {params.unzipped} {params.interm1} > {params.interm2}
#            awk 'BEGIN{{OFS="	"}} {{print $13, $14, $1, $7, $9, $11, $12}}' {params.interm2} > {output}
#            rm {params.unzipped} {params.interm1} {params.interm2}
#        '''
#
#
##rules.parse_mu_et_al_eqtl.output:
#rule run_colocalization_mu_et_al_example:
#    input:
#        gwas = 'results/main/coloc/Data/T1D_GWAS/T1D_34012112_Gaulton/GWAS_input_colocalization_pval_lt_5eMinus8.txt',
#        eqtl = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/BLUEPRINT_eQTL/Monocyte.txt.gz'
#    output:
#        outdir = directory('results/main/coloc/Results/Colocalization_SMK/T1D_34012112_Gaulton/BLUEPRINT_eQTL_Monocyte/')
#    params:
#        chr = 1,
#        pos = 2,
#        gene_name = 3,
#        dist = 4,
#        slope = 5,
#        pvalue = 6,
#        fdr = 7,
#    resources:
#        mem_mb = 20000,
#        nodes = 1,
#        ppn = 1,
#        #PBS -m ae
#        #PBS -j eo
#        #PBS -V
#    shell:
#        r"""
#            TMPDIR=/scratch
#            #cd {output}
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
#                                --eqtl-header \
#                                {input.gwas} \
#                                {input.eqtl} \
#                                {output.outdir} 
#        """

# This file doesn not contain a chr and position field, rather it has a 
# sid filed form which you can generated these two.
rule add_dist_col:
    input:
        gencode = 'results/refs/gencode/v39/gencode.v39lift37.annotation.bed',
        eqtl = 'results/main/eqtl/{eqtl_source}/ge/{eqtl_source}_ge_{ge_source}.all.tsv.gz'
    output:
        'results/main/eqtl/{eqtl_source}/ge/{eqtl_source}_ge_{ge_source}.all.dist.tsv.gz'
    log: 
        'results/main/eqtl/logs/add_dist_col.{eqtl_source}.{ge_source}.log'
    resources:
        mem_mb = 16000,
        nodes = 1,
        ppn = 1,
    shell:
        r'''
            /mnt/BioHome/jreyna/software/anaconda3/envs/hic_tls/bin/python \
                    workflow/scripts/eqtl/Add_Dist_Col.py \
                        --header \
                        {input.gencode} \
                        {input.eqtl} \
                        {output} >> {log} 2>&1
        '''

# This file doesn not contain a chr and position field, rather it has a 
# sid filed form which you can generated these two.
rule parse_eqtl_catalog_with_dist:
    input:
        rules.add_dist_col.output
    output:
        'results/main/eqtl/{eqtl_source}/ge/{eqtl_source}_ge_{ge_source}.all.dist.input.tsv.gz'
    params:
        interm='results/main/eqtl/{eqtl_source}/ge/{eqtl_source}_ge_{ge_source}.all.dist.input.tsv'
    log: 
        'results/main/coloc/Results/Colocalization_SMK/logs/parse_eqtl_catalog_with_dist.{eqtl_source}.{ge_source}.log'
    shell:
        r'''
            zcat {input} | \
                 sed '1d' | \
                 awk 'BEGIN{{OFS="	"}}; {{print $2, $3, $17, $20, $10, $9}}' > \
                 {params.interm} 2> {log}
            bgzip {params.interm}
        '''
    











# This file doesn not contain a chr and position field, rather it has a 
# sid filed form which you can generated these two.
rule parse_immunexut_eqtl:
    input:
        'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.txt.gz'
    output:
        'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.input.txt'
    params:
        unzipped = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.txt',
        interm1 = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.interm1.txt',
        interm2 = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.interm2.txt'
    log: 
        'results/main/coloc/Results/Colocalization_SMK/logs/parse_mu_et_al_eqtl.{eqtl_source}.{ge_source}.log'
    shell:
        r'''
            zcat {input} > {params.unzipped}
             cut -f 6 {params.unzipped} | \
                 sed '1d' | \
                 awk 'BEGIN{{FS=":"; OFS=" "; print "chr   pos"}} {{print "chr"$1,$2}}' > {params.interm1}
            paste {params.unzipped} {params.interm1} > {params.interm2}
            awk 'BEGIN{{OFS="	"}} {{print $13, $14, $1, $7, $9, $11, $12}}' {params.interm2} > {output}
            rm {params.unzipped} {params.interm1} {params.interm2}
        '''






# This file doesn not contain a chr and position field, rather it has a 
# sid filed form which you can generated these two.
rule parse_mu_et_al_eqtl:
    input:
        'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.txt.gz'
    output:
        'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.input.txt'
    params:
        unzipped = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.txt',
        interm1 = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.interm1.txt',
        interm2 = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.interm2.txt'
    log: 
        'results/main/coloc/Results/Colocalization_SMK/logs/parse_mu_et_al_eqtl.{eqtl_source}.{ge_source}.log'
    shell:
        r'''
            zcat {input} > {params.unzipped}
             cut -f 6 {params.unzipped} | \
                 sed '1d' | \
                 awk 'BEGIN{{FS=":"; OFS=" "; print "chr   pos"}} {{print "chr"$1,$2}}' > {params.interm1}
            paste {params.unzipped} {params.interm1} > {params.interm2}
            awk 'BEGIN{{OFS="	"}} {{print $13, $14, $1, $7, $9, $11, $12}}' {params.interm2} > {output}
            rm {params.unzipped} {params.interm1} {params.interm2}
        '''

#rules.parse_mu_et_al_eqtl.output:
rule run_colocalization_mu_et_al:
    input:
        gwas = 'results/main/coloc/Data/T1D_GWAS/{gwas_source}/GWAS_input_colocalization_pval_lt_5eMinus8.txt',
        eqtl = 'results/main/coloc/Data/eqtl_sqtl_summ_stats/{eqtl_source}/{ge_source}.input.txt'
    output:
        outdir = directory('results/main/coloc/Results/Colocalization_SMK/{gwas_source}/{eqtl_source}/{ge_source}/')
    log: 
        'results/main/coloc/Results/Colocalization_SMK/logs/run_colocalization_mu_et_al_example.{gwas_source}.{eqtl_source}.{ge_source}.log'
    params:
        chr = 1,
        pos = 2,
        gene_name = 3,
        dist = 4,
        slope = 5,
        pvalue = 6,
        fdr = 7,
        header = 'TRUE'
    resources:
        mem_mb = 20000,
        nodes = 1,
        ppn = 1,
        #PBS -j eo
        #PBS -V
    shell:
        r"""
            # make the log directory
            mkdir -p $(dirname {log})

            # set the tempdir 
            TMPDIR=/scratch

            # add the corresponding R for colocalization analyses
            PATH=/mnt/BioApps/R/3.6.1/bin/:$PATH

            # run colocalization command
            Rscript scripts/coloc/Colocalization_Analysis_GWAS_Script_Generalized.R \
                                --eqtl-chr {params.chr} \
                                --eqtl-pos {params.pos} \
                                --eqtl-geneName {params.gene_name} \
                                --eqtl-dist {params.dist} \
                                --eqtl-slope {params.slope} \
                                --eqtl-pvalue {params.pvalue} \
                                --eqtl-FDR {params.fdr} \
                                --eqtl-header {params.header} \
                                {input.gwas} \
                                {input.eqtl} \
                                {output.outdir} >> {log} 2>&1
        """
    

