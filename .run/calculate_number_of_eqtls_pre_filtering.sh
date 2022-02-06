#########################################################################################
# Getting all eQTL datasets
#########################################################################################

outfiles+=""
for fn in $(ls results/main/eqtl/*/ge/*_ge_*.all.tsv.gz);
do 
    # extracting the eqtl_source and ge_source
    keynames=$(echo $fn | xargs basename | cut -d "_" -f 1,3 --output-delimiter "	")
    eqtl_source=$(echo $fn | xargs basename | sed "s/_ge_/	/" | cut -d "	" -f 1)
    ge_source=$(echo $fn | xargs basename | sed "s/_ge_/	/" | cut -d "	" -f 2 | cut -f 1 -d .)

    #echo "keynames: $keynames"
    echo "eqtl_source: $eqtl_source"
    echo "ge_source: $ge_source"

    # getting the output files
    new_fn="results/main/eqtl/${eqtl_source}/ge/${eqtl_source}_ge_${ge_source}.all.prefilter.num_eqtls.txt"
    outfiles+=" $new_fn"
    echo "$fn ----------> $new_fn"
done

#outfiles=$(echo $outfiles | cut -d " " -f 1)
echo
echo "####################################################################################"
echo "outfiles: $outfiles"
echo "####################################################################################"

##snakemake --profile workflow/profiles/pbs-torque/ -n $outfiles
snakemake --profile workflow/profiles/pbs-torque/ $@ $outfiles
#snakemake --profile workflow/profiles/local/ $@ $outfiles
#echo "snakemake --profile workflow/profiles/local/ $@ $outfiles"
