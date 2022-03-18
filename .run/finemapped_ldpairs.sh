fns=""
while IFS= read -r line; 
do
    run_info=($line)
    fn="results/main/ldpairs/finemapping/${run_info[0]}/GRCh37/offset_1000000/finemapping_ld_snps.txt"
    fns+="$fn "
done < <(cat config/gwas.samplesheet.tsv | sed '1d' | grep -v "^#")
#snakemake --profile workflow/profiles/local $@ $fns
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
