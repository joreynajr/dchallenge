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

        

rule filter_eqtls_using_colocs: #(Status: developing)
    input:
        eqtl = rules.add_missing_cols.output,
        coloc_dir = rules.run_colocalization_eqtl_catalog.output,
        coloc = 'results/main/coloc/Results/Colocalization_SMKN/T1D_34012112_Gaulton/{eqtl_source}/{ge_source}/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed'
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


rule annotate_colocs: #(Status: developing)
    input:
        eqtl = rules.filter_eqtls_using_colocs.output,
        coloc_dir = rules.run_colocalization_eqtl_catalog.output,
        coloc = 'results/main/coloc/Results/Colocalization_SMKN/T1D_34012112_Gaulton/{eqtl_source}/{ge_source}/FINAL_Summary_Coloc_Gene_SNP_Pairs.bed',
        loops = 'results/main/h3k27ac_hichip/{ge_source}/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.bed.gz',
        genome_sizes = 'results/refs/hg38/hg38.chrom.sizes',
        gencode = 'results/refs/gencode/v30/gencode.v30.annotation.bed'
    output:
        outdir = directory('results/main/sgls/{gwas_source}/{eqtl_source}/{ge_source}/{loop_source}/')
    log: 
        'results/main/sgls/logs/annotate_colocs.{gwas_source}.{eqtl_source}.{ge_source}.{loop_source}.log'
    shell:
        r'''
            mkdir -p {output}
            python workflow/scripts/sgls/Find_SGLs.eQTL_Catalogue_Format.py \
                                        {input.eqtl} \
                                        {input.coloc} \
                                        {input.loops} \
                                        {wildcards.ge_source} \
                                        {input.genome_sizes} \
                                        {input.gencode} \
                                        {output}
        '''
