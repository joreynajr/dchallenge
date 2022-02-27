rule run_finemapping:
    input:
        'results/main/finemapping/{gwas_source}/{ref_genome}/raw/GWAS_RawData.finemap.input.tsv'
    output:
        directory('results/main/finemapping/{gwas_source}/{ref_genome}/offset_{offset}/')
    params:
        chrcol = 1,
        poscol = 2,
        pvalcol = 3,
        afcol = 4,
        betacol = 5,
        secol = 6 
    resources:
        mem_mb = 10000,
        nodes = 1,
        ppn = 8
    shell:
        r'''
            # make a config file
            mkdir -p {output}
            abs_config=$(readlink -f "{output}/config.txt")
            abs_input=$(readlink -f {input})
            abs_output=$(readlink -f {output})
            cat workflow/scripts/finemapping/Finemap_SMK/configfile.template.txt | \
                        sed "s+{{gwas}}+$abs_input+" | \
                        sed "s+{{outdir}}+$abs_output+" | \
                        sed "s+{{offset}}+{wildcards.offset}+" | \
                        sed "s+{{chrcol}}+{params.chrcol}+" | \
                        sed "s+{{poscol}}+{params.poscol}+" | \
                        sed "s+{{pvalcol}}+{params.pvalcol}+" | \
                        sed "s+{{afcol}}+{params.afcol}+" | \
                        sed "s+{{betacol}}+{params.betacol}+" | \
                        sed "s+{{secol}}+{params.secol}+"  > $abs_config

            # run the finemapping
            cd workflow/scripts/finemapping/Finemap_GWAS_Joaquin/
            bash Finemap_Script.sh -C $abs_config
        '''
