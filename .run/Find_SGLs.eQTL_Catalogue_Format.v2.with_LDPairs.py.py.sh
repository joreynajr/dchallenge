python="/mnt/BioAdHoc/Groups/vd-ay/jreyna/software/mambaforge/envs/hichip-db/bin/python"
loop_slop=25000
while IFS= read -r line;
do
    # reading in parameters
    run_info=($line)

    #echo "run_info ${run_info[@]}"

    gwas_source="${run_info[0]}"
    ref_genome="GRCh37"
    eqtl_source="${run_info[1]}"
    ge_source="${run_info[2]}"
    loop_source="${run_info[3]}"

    # setting paths 
    eqtl="results/main/sgls/${gwas_source}/${eqtl_source}/${ge_source}/eqtls.coloc_filtered.tsv.gz"
    coloc="results/main/coloc/eQTL_Catalogue/${gwas_source}/${ref_genome}/${eqtl_source}/${ge_source}/ldpairs/coloc_ld_snps.txt"
    loops="results/main/h3k27ac_hichip/${loop_source}/FitHiChIP_S/FitHiChIP.interactions_FitHiC_Q0.01.grch38.bed"
    genome_sizes="results/refs/hg38/hg38.chrom.sizes"
    gencode="results/refs/gencode/v30/gencode.v30.annotation.bed"
    outdir="results/main/sgls/ldpairs/${gwas_source}/${ref_genome}/${eqtl_source}/${ge_source}/${loop_source}/"
    log="results/main/sgls/logs/annotate_colocs.${gwas_source}.${ref_genome}.${eqtl_source}.${ge_source}.${loop_source}.log"
    celltype="$ge_source"

    # make the outut file
    mkdir -p $outdir 

    #if [ ! -d "$outdir" ];
    #then
        # #valute the current rule
        #cmd="$python workflow/scripts/sgls/Find_SGLs.eQTL_Catalogue_Format.v2.with_LDPairs.py.py"
        #cmd+=" $eqtl $coloc $loops $celltype $genome_sizes $gencode $loop_slop $outdir"
        #echo "cmd: $cmd"
        #eval $cmd | tee $log
    #fi

    # evalute the current rule
    cmd="$python workflow/scripts/sgls/Find_SGLs.eQTL_Catalogue_Format.v2.with_LDPairs.py.py"
    cmd+=" $eqtl $coloc $loops $celltype $genome_sizes $gencode $loop_slop $outdir"
    echo "cmd: $cmd"
    eval $cmd | tee $log
    
done < <(cat config/sgl_samplesheets/sgl.samplesheet.t1d_working.tsv | sed '1d' | grep -v "^#" | head -n 1000)
