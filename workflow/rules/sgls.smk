# Rename the loop directories for easier access
#'results/main/2021_Nikhil_eQTL/Data/FitHiChIP_Loops/CD4N/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed'
rule rename_loop_dirs: #(Status: running)
    input:
        indir = 'results/main/2021_Nikhil_eQTL/Data/FitHiChIP_Loops/',
        fithic2eqtl_cat = 'results/refs/fithic_to_eqtl_cat_names.tsv'
    output:
        outdir = directory('results/main/h3k27ac_hichip/')
    shell:
        r'''
            mkdir -p {output}
            old_names=($(cut -f 1 {input.fithic2eqtl_cat}))
            new_names=($(cut -f 2 {input.fithic2eqtl_cat}))
            for i in ${{!old_names[@]}};
            do
                # get the names
                old_name=${{old_names[$i]}}
                new_name=${{new_names[$i]}}

                # get the directory names
                old_dir="{input.indir}/${{old_name}}"
                new_dir="{output}/${{new_name}}"

                # print the directory names
                echo "old_dir: $old_dir"
                echo "new_dir: $new_dir"

                # if the old_dir is present then move it
                if [[ -d $old_dir ]];
                then
                    ln -sr $old_dir $new_dir
                fi
            done
        '''


# Filter for eQTLs with colocalizations for annotations purposes
# Issue: Doesnt seem like I need it because the colocalization script takes 
# care of this now. 
rule filter_eqtls_using_colocs: #(Status: developing)
    input:
        eqtl = rules.add_missing_cols.output,
        coloc_dir = rules.run_colocalization_eqtl_catalog.output,
        coloc = 'results/main/coloc/Results/eQTL_Catalogue/{gwas_source}/{eqtl_source}/{ge_source}/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed'
    output:
        outfn = 'results/main/sgls/{gwas_source}/{eqtl_source}/{ge_source}/eqtls.coloc_filtered.tsv.gz'
    log: 
        'results/main/sgls/logs/filter_eqtls_using_colocs.{gwas_source}.{eqtl_source}.{ge_source}.log'
    shell:
        r'''
            python workflow/scripts/sgls/Filter_eQTLs_with_Colocs.py \
                                        {input.eqtl} \
                                        {input.coloc} \
                                        {output}
        '''


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


rule combine_intersections: #(Status: developing)
    output:
        'results/main/sgls/combined/super_master.snp_gene_loop.analysis.tsv'
    log: 
        'results/main/sgls/logs/combine_intersections.log'
    notebook:
        '../notebooks/reports/Combine_Master_Tables.py.ipynb'
