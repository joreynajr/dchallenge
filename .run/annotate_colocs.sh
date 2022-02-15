fns=""
while IFS= read -r line;
do
    run_info=($line)

    #echo "${run_info[0]}/${run_info[1]}/${run_info[2]}/${run_info[3]}/"
    fn="results/main/sgls/${run_info[0]}/${run_info[1]}/${run_info[2]}/${run_info[3]}/"
    dep_fn="results/main/coloc/Results/Colocalization_SMKN/${run_info[0]}/${run_info[1]}/${run_info[2]}"
    if [[ -d "$dep_fn" ]];
    then
        echo "PROCESSING: $fn"
        fns+="$fn "
    else
        echo "WARNING: Missing dependencies for $fn"
    fi
done < <(cat config/sgl.samplesheet.tsv.v2 | sed '1d' | grep -v "^#" | head -n 100)
#done < <(cat config/sgl.samplesheet.tsv | sed '1d')

echo "fns: $fns"
#snakemake --profile workflow/profiles/local/ $@ $fns
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns

