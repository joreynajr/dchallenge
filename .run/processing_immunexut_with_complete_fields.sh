#########################################################################################
# Getting all eQTL datasets
#########################################################################################

eqtl_source="ImmuNexUT"
outfiles+=""
for fn in $(ls results/main/eqtl/ImmuNexUT/raw/E-GEAD-420/*_nominal.txt);
do 
    ge_source=$(basename $fn | sed 's/_nominal\.txt//')

    # getting the output files
    new_fn="results/main/eqtl/${eqtl_source}/ge/${eqtl_source}_ge_${ge_source}.all.immunexut.complete_fields.input.tsv.gz"
    outfiles+=" $new_fn"
    echo "$fn ----------> $new_fn"
done

#echo $outfiles

#outfiles="results/main/eqtl/BLUEPRINT/ge/BLUEPRINT_ge_monocyte.all.complete_fields.input.tsv.gz"

##snakemake --profile workflow/profiles/pbs-torque/ -n $outfiles
snakemake --profile workflow/profiles/pbs-torque/ $@ $outfiles
#snakemake --profile workflow/profiles/local/ $@ $outfiles
#echo "snakemake --profile workflow/profiles/local/ $@ $outfiles"
