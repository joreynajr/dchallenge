
fns=""
ref_genome="GRCh37"
offset=1000000
for gwas in $(cut -f 1 config/gwas.samplesheet.tsv | sed '1d' | grep -v "^#" | head -n 1);
do
    fn="results/main/finemapping/${gwas}/${ref_genome}/offset_${offset}/"
    fns+="$fn "
done

#echo $fns | sed 's/ /\n/g'
#snakemake --profile workflow/profiles/local $@ $fns
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
