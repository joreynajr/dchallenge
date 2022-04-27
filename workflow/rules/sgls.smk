# Filter for eQTLs with colocalizations for annotations purposes
# Issue: Doesnt seem like I need it because the colocalization script takes 
# care of this now. 
#
def get_filter_eqtls_using_colocs_eqtl(wildcards): #(Status: running)
    if wildcards.eqtl_source == 'ImmuNexUT':
        eqtl_fn = 'results/main/eqtl/ImmuNexUT/ge/ImmuNexUT_ge_{ge_source}.all.immunexut.dist.fdr.tsv.gz'
    else:
        eqtl_fn = 'results/main/eqtl/{eqtl_source}/ge/{eqtl_source}_ge_{ge_source}.all.dist.fdr.tsv.gz'
    eqtl_fn = eqtl_fn.format(**wildcards)
    print('eqtl_fn: {}'.format(eqtl_fn))
    return(eqtl_fn)


# eqtl = rules.add_missing_cols.output,
rule filter_eqtls_using_colocs: #(Status: running)
    input:
        eqtl = get_filter_eqtls_using_colocs_eqtl,
        coloc = 'results/main/coloc/Results/{eqtl_db}/{gwas_source}/{eqtl_source}/{ge_source}/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed'
    output:
        outfn = 'results/main/sgls/{eqtl_db}/{gwas_source}/{eqtl_source}/{ge_source}/eqtls.coloc_filtered.tsv.gz'
    log: 
        'results/main/sgls/logs/filter_eqtls_using_colocs.{eqtl_db}.{gwas_source}.{eqtl_source}.{ge_source}.log'
    shell:
        r'''
            python workflow/scripts/sgls/Filter_eQTLs_with_Colocs.py \
                                        {input.eqtl} \
                                        {input.coloc} \
                                        {output}
        '''


# liftover the significant eQTLs to GRCh37
rule liftover_sig_eqtls_to_GRCh37: # (Status: running)
    input:
        rules.filter_eqtls_using_colocs.output.outfn
    output:
        outfn = 'results/main/GRCh37/sgls/{eqtl_db}/{gwas_source}/{eqtl_source}/{ge_source}/eqtls.coloc_filtered.tsv.gz'
    params:
        chr = 2,
        pos = 3,
        sep = '\t',
        header = 'TRUE',
        temp_gunzip = 'results/main/GRCh37/sgls/{eqtl_db}/{gwas_source}/{eqtl_source}/{ge_source}/eqtls.coloc_filtered.temp.tsv',
        temp_liftover = 'results/main/GRCh37/sgls/{eqtl_db}/{gwas_source}/{eqtl_source}/{ge_source}/eqtls.coloc_filtered.tsv'
    log: 
        'results/main/GRCh37/sgls/logs/liftover_sig_eqtls_to_GRCh37.{eqtl_db}.{gwas_source}.{eqtl_source}.{ge_source}.log'
    shell:
        r"""
            # uncompress the input
            gzip -d -c {input} > {params.temp_gunzip} 2> {log}

            # liftover from GRCh38 to GRCh37
            {config[hichip_db_py]} workflow/scripts/utilities/general_liftover.py \
                    -i {params.temp_gunzip} \
                    -o {params.temp_liftover} \
                    --chr-col {params.chr} \
                    --pos-col {params.pos} \
                    --header \
                    --sep "{params.sep}" > {log} 2>&1

            # compress the final temp file
            gzip {params.temp_liftover} > {log} 2>&1

        """


# Liftover loops from GRCh37 to GRCh38
# Needs and update to use on the start
rule liftover_loops_to_GRCh38: #(Status: running)
    input:
        loops = 'results/main/h3k27ac_hichip/{loop_source}/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed',
        chain = rules.download_chain_file_hg19tohg38.output
    output:
        loops = 'results/main/h3k27ac_hichip/{loop_source}/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.grch38.bed'
    log: 
        'results/main/sgls/logs/liftover_loops_to_GRCh38.{loop_source}.log'
    resources:
        mem_mb = 24000,
        nodes = 1,
        ppn = 1,
    shell:
        r'''
            # lift the left anchor
            echo "# lift the left anchor" > {log}
            left="{input.loops}.left"
            awk 'BEGIN{{OFS="	"}} {{if(NR != 1) {{print $1, $2, $3, NR}}}}' {input.loops} > $left 2> {log}
            left_lifted="{input.loops}.left.lifted"
            left_unmapped="{input.loops}.left.unmapped"
            {config[liftover]} -bedPlus=3 -tab $left {input.chain} $left_lifted $left_unmapped 2> {log}

            # sort the left data
            echo "# sort the left data" > {log}
            left_sorted="{input.loops}.left.lifted.sorted"
            sort -k 4 $left_lifted > $left_sorted 2> {log}

            # lift the right anchor
            echo "# lift the right anchor" > {log}
            right="{input.loops}.right"
            sed '1d' {input.loops} | cut -f 4- | awk 'BEGIN{{OFS="	"}} {{print $0, NR}}' > $right 2> {log}
            right_lifted="{input.loops}.right.lifted"
            right_unmapped="{input.loops}.right.unmapped"
            {config[liftover]} -bedPlus=3 -tab $right {input.chain} $right_lifted $right_unmapped 2> {log}

            # sort the right data
            echo "# sort the right data" > {log}
            right_id=$(head -n 1 $right_lifted | wc -w)
            right_sorted="{input.loops}.right.lifted.sorted"
            sort -k $right_id $right_lifted > $right_sorted 2> {log}

            # join the two files on the serial ID
            head -n 1 {input.loops} > {output} 2> {log} # add header fist
            join -1 4 -2 $right_id -t "	" $left_sorted $right_sorted | cut -f 2- >> {output} 2> {log}

            # remove temps
            rm $left_lifted $left_unmapped $left_sorted 2> {log}
            rm $right_lifted $right_unmapped $right_sorted 2> {log}

        '''


## Annotation for colocalized SNP-Gene pairs using loops
#rule annotate_colocs: #(Status: developing; in the crosshairs of deprecation)
#    input:
#        eqtl = rules.filter_eqtls_using_colocs.output,
#        coloc_dir = rules.run_colocalization_eqtl_catalog.output,
#        coloc = 'results/main/coloc/Results/Colocalization_SMKN/{gwas_source}/{eqtl_source}/{ge_source}/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed',
#        loops = rules.liftover_loops_to_GRCh38.output.loops,
#        genome_sizes = 'results/refs/hg38/hg38.chrom.sizes',
#        gencode = 'results/refs/gencode/v30/gencode.v30.annotation.bed'
#    output:
#        outdir = directory('results/main/sgls/{gwas_source}/{eqtl_source}/{ge_source}/{loop_source}/')
#    log: 
#        'results/main/sgls/logs/annotate_colocs.{gwas_source}.{eqtl_source}.{ge_source}.{loop_source}.log'
#    params:
#        loop_slop = 25000
#    resources:
#        mem_mb = 24000,
#        nodes = 1,
#        ppn = 1,
#    shell:
#        r'''
#            mkdir -p {output} >> {log} 2>&1
#            {config[hichip_db_py]} workflow/scripts/sgls/Find_SGLs.eQTL_Catalogue_Format.py \
#                                        {input.eqtl} \
#                                        {input.coloc} \
#                                        {input.loops} \
#                                        {wildcards.ge_source} \
#                                        {input.genome_sizes} \
#                                        {input.gencode} \
#                                        {params.loop_slop} \
#                                        {output} >> {log} 2>&1
#        '''


## Annotation for colocalized SNP-Gene pairs using loops
## and producing a notebook as a log file
#rule annotate_colocs_notebook: #(Status: developing)
#    input:
#        eqtl = rules.filter_eqtls_using_colocs.output,
#        coloc_dir = rules.run_colocalization_eqtl_catalog.output,
#        coloc = 'results/main/coloc/Results/eQTL_Catalogue/{gwas_source}/{eqtl_source}/{ge_source}/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed',
#        loops = rules.liftover_loops_to_GRCh38.output.loops,
#        genome_sizes = 'results/refs/hg38/hg38.chrom.sizes',
#        gencode = 'results/refs/gencode/v30/gencode.v30.annotation.bed'
#    output:
#        outdir = directory('results/main/sgls/{gwas_source}/{eqtl_source}/{ge_source}/{loop_source}/')
#    log: 
#        #notebook='results/main/sgls/logs/annotate_colocs.{gwas_source}.{eqtl_source}.{ge_source}.{loop_source}.ipynb' #not working, JSON Error
#        'results/main/sgls/logs/annotate_colocs.{gwas_source}.{eqtl_source}.{ge_source}.{loop_source}.log'
#    params:
#        loop_slop = 25000
#    resources:
#        mem_mb = 24000,
#        nodes = 1,
#        ppn = 1,
#    notebook:
#            "../scripts/sgls/Find_SGLs.eQTL_Catalogue_Format.v2.py.ipynb"



# Annotation for colocalized SNP-Gene pairs using loops
# and producing a notebook as a log file
rule annotate_colocs_script: #(Status: developing)
    input:
        eqtl = rules.filter_eqtls_using_colocs.output,
        coloc_dir = rules.run_colocalization_eqtl_catalog.output,
        coloc = 'results/main/coloc/Results/eQTL_Catalogue/{gwas_source}/{eqtl_source}/{ge_source}/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed',
        loops = rules.liftover_loops_to_GRCh38.output.loops,
        genome_sizes = 'results/refs/hg38/hg38.chrom.sizes',
        gencode = 'results/refs/gencode/v30/gencode.v30.annotation.bed'
    output:
        outdir = directory('results/main/sgls/{gwas_source}/{eqtl_source}/{ge_source}/{loop_source}/')
    log: 
        #notebook='results/main/sgls/logs/annotate_colocs.{gwas_source}.{eqtl_source}.{ge_source}.{loop_source}.ipynb' #not working, JSON Error
        'results/main/sgls/logs/annotate_colocs.{gwas_source}.{eqtl_source}.{ge_source}.{loop_source}.log'
    params:
        loop_slop = 25000
    resources:
        mem_mb = 24000,
        nodes = 1,
        ppn = 1,
    script:
            "../scripts/sgls/Find_SGLs.eQTL_Catalogue_Format.v2.py"

# combine all annotations from rule annotate_colocs_script
rule combine_intersections: #(Status: developing)
    output:
        'results/main/sgls/combined/super_master.snp_gene_loop.analysis.tsv'
    log: 
        'results/main/sgls/logs/combine_intersections.log'
    notebook:
        '../notebooks/reports/Combine_Master_Tables.py.ipynb'


# DEPRECATED
## Annotation for colocalized SNP-Gene pairs using loops
## and producing a notebook as a log file
##loops = rules.liftover_loops_to_GRCh38.output.loops,
#rule Find_SGL_for_Coloc_and_LD: #(Status: developing)
#    input:
#        eqtl = rules.liftover_sig_eqtls_to_GRCh37.output.outfn,
#        coloc = rules.find_ldpairs_for_coloc_snps.output.ld,
#        loops = 'results/main/h3k27ac_hichip/{loop_source}/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed',
#        genome_sizes = 'results/refs/hg19/hg19.chrom.sizes',
#        gencode = 'results/refs/gencode/v30/gencode.v30.annotation.grch37.bed'
#    output:
#        outdir = directory('results/main/GRCh37/sgls/ldpairs/{gwas_source}/{eqtl_source}/{ge_source}/{loop_source}/')
#    log: 
#        #notebook='results/main/sgls/logs/annotate_colocs.{gwas_source}.{eqtl_source}.{ge_source}.{loop_source}.ipynb' #not working, JSON Error
#        'results/main/sgls/logs/annotate_colocs.{gwas_source}.{eqtl_source}.{ge_source}.{loop_source}.log'
#    params:
#        loop_slop = 25000
#    resources:
#        mem_mb = 24000,
#        nodes = 1,
#        ppn = 1,
#    notebook:
#            "../scripts/sgls/Find_SGLs.eQTL_Catalogue_Format.v2.with_LDPairs.py.ipynb"

# Annotation for colocalized SNP-Gene pairs using loops
# and producing a notebook as a log file
#loops = rules.liftover_loops_to_GRCh38.output.loops,
rule Find_SGL_for_Coloc_and_LD_Script_Version: #(Status: developing)
    input:
        eqtl = rules.liftover_sig_eqtls_to_GRCh37.output.outfn,
        coloc = rules.find_ldpairs_for_coloc_snps.output.ld,
        loops = 'results/main/h3k27ac_hichip/{loop_source}/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed',
        genome_sizes = 'results/refs/hg19/hg19.chrom.sizes',
        gencode = 'results/refs/gencode/v30/gencode.v30.annotation.grch37.bed'
    output:
        outdir = directory('results/main/GRCh37/sgls/ldpairs/{eqtl_db}/{gwas_source}/{eqtl_source}/{ge_source}/{loop_source}/script_version/')
    log: 
        #notebook='results/main/sgls/logs/annotate_colocs.{gwas_source}.{eqtl_source}.{ge_source}.{loop_source}.ipynb' #not working, JSON Error
        'results/main/sgls/logs/annotate_colocs.{eqtl_db}.{gwas_source}.{eqtl_source}.{ge_source}.{loop_source}.log'
    params:
        loop_slop = 25000
    resources:
        mem_mb = 24000,
        nodes = 1,
        ppn = 1,
    shell:
        r"""
            {config[hichip_db_py]} workflow/scripts/sgls/Find_SGLs.eQTL_Catalogue_Format.v3.with_LDPairs.py \
                {input.eqtl} \
                {input.coloc} \
                {input.loops} \
                {wildcards.eqtl_source} \
                {input.genome_sizes} \
                {input.gencode} \
                {params.loop_slop} \
                {output} > {log} 2>&1
        """
