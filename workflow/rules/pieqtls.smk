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
