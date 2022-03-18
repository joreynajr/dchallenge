#########################################################################################
# Getting all eQTL datasets
#########################################################################################

# Replace by new approach
#outfiles+=""
#for fn in $(ls results/main/eqtl/*/ge/*_ge_*.all.tsv.gz);
#do 
#    # extracting the eqtl_source and ge_source
#    keynames=$(echo $fn | xargs basename | cut -d "_" -f 1,3 --output-delimiter "	")
#    eqtl_source=$(echo $keynames | cut -d " " -f 1)
#    ge_source=$(echo $keynames | cut -d " " -f 2 | cut -d . -f 1)
#
#    echo $eqtl_source $ge_source
#
#    #echo "keynames: $keynames"
#    #echo "eqtl_source: $eqtl_source"
#    #echo "ge_source: $ge_source"
#
#    # getting the output files
#    new_fn="results/main/eqtl/${eqtl_source}/ge/${eqtl_source}_ge_${ge_source}.all.complete_fields.input.tsv.gz"
#    outfiles+=" $new_fn"
#    echo "$fn ----------> $new_fn"
#done

fns=""
while IFS= read -r line;
do
    run_info=($line)
    fn="results/main/eqtl/${run_info[0]}/ge/${run_info[0]}_ge_${run_info[1]}.all.dist.fdr.input.tsv.gz"
    fns+="$fn "

    #dep_fn="results/main/coloc/Results/eQTL_Catalogue/${run_info[0]}/${run_info[1]}/${run_info[2]}"
    #if [[ -d "$dep_fn" ]];
    #then
    #    echo "PROCESSING: $fn"
    #    fns+="$fn "
    #else
    #    echo "WARNING: Missing dependencies $dep_fn for $fn"
    #fi
done < <(cat config/sgl.samplesheet.tsv.v2 | sed '1d' | grep -v "^#" | cut -f 2,3 | sort | uniq | head -n 100)

##snakemake --profile workflow/profiles/pbs-torque/ -n $outfiles
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
#snakemake --profile workflow/profiles/local/ $@ $outfiles
#echo "snakemake --profile workflow/profiles/local/ $@ $outfiles"
