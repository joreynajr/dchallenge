fns=""
while IFS= read -r line;
do
    run_info=($line)
    fns+="results/main/sgls/${run_info[0]}/${run_info[1]}/${run_info[2]}/${run_info[3]}/ "
    echo $fns
    break
done < <(cat config/sgl.samplesheet.tsv | sed '1d')

echo "fns: $fns"
snakemake --profile workflow/profiles/local/ $@ $fns

