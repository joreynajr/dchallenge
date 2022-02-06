# These are general purpose rules that are commonly reused across the pipeline

# convert from fithic_to_bedpe
rule fithic_to_bedpe:
    shell:
        r"""
           {config[hic_tls_py]} workflow/scripts/fithic_format_to_bedpe.py \
                            --fithic {input} \
                            --bedpe {output} \
                            --res {wildcards.res} \
                            --qmax {params.qmax} >> {log} 2>&1
        """

# annotate bedpe for use with JuiceBox 
rule bedpe_to_juicer_2d_annotations:
    conda:
        '../envs/environment.yml'
    shell:
        """
            {config[hic_tls_py]} workflow/scripts/bedpe_to_annotated_bedpe.py \
                    --bedpe {input.bedpe} \
                    --anno {output.anno} \
                    --color {params.color} >> {log} 2>&1
        """

# Generating the input data to create a HIC file from a validPairs file
# Must define the following:
#input:
#validpair = 'data/hic_validpairs/{cline}.{srr_id}_GCA_000001405.15.GRCh38_no_alt_analysis_set.fna.bowtie_index.bwt2pairs.validPairs',
#    genome_sizes = rules.download_hg38_files.output.genome_sizes,
#params:
#    jt_input = 'results//main/{cline}/juicebox/{cline}.{srr_id}.juicertools.input.txt'
#output:
#    hic = 'results/main/{cline}/preprocessing/{cline}.{srr_id}.hic'
#log:
#    'results/main/{cline}/logs/rule_validpairs_to_juicertools_input_raw_matrix_counts_{cline}_{srr_id}.out'
#resources: # has been set here but can be modified in child rule
#    mem_mb = 80000
rule validpairs_to_hic:
    resources:
        mem_mb = 80000
    shell:
        """
            # convert validpairs to juicertools input
            awk '{{
                    if($4=="+") {{s1=0}} else {{s1=1}}; \
                    if($7=="+") {{s2=0}} else {{s2=1}}; \
                    if($2 < $5) {{print $2" "$5" "s1" "s2" "$3" "$6" 0 1 "$1}} \
                    else{{print $5" "$2" "s2" "s1" "$6" "$3" 0 1 "$1}}}}' {input} | \
           sort -k 1,1 -k2,2n | \
           awk '{{print $9"\t"$3"\t"$1"\t"$5"\t"$7"\t"$4"\t"$2"\t"$6"\t"$8"\t42\t42"}}' > {params.jt_input} 2> {log}

            # convert juicertools input to hic
            java -jar {config[juicer_tools]} pre {params.jt_input} {output.hic} {input.genome_sizes} >> {log} 2>&1
        """
