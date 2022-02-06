#########################################################################################
# Getting all eQTL datasets
#########################################################################################

eqtl_source="ImmuNexUT"
outfiles+=""
for fn in $(ls results/main/eqtl/ImmuNexUT/raw/E-GEAD-420/*_nominal.txt);
do 
    ge_source=$(basename $fn | sed 's/_nominal\.txt//')
    new_fn="results/main/eqtl/${eqtl_source}/ge/${eqtl_source}_ge_${ge_source}.all.immunexut.fdr.txt"
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
