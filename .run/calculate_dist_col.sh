#########################################################################################
# Getting all eQTL datasets
#########################################################################################

outfiles+=""
for fn in $(ls results/main/eqtl/*/ge/*_ge_*.all.tsv.gz);
do 
    new_fn=$(echo $fn | sed 's/\(\.tsv\.gz\)/.dist.txt/')
    outfiles+=" $new_fn"
    echo "$fn ----------> $new_fn"
done

#outfiles=$(echo $outfiles | cut -d " " -f 2)
echo
echo "####################################################################################"
echo "outfiles: $outfiles"
echo "####################################################################################"

##snakemake --profile workflow/profiles/pbs-torque/ -n $outfiles
snakemake --profile workflow/profiles/pbs-torque/ $@ $outfiles
#snakemake --profile workflow/profiles/local/ $@ $outfiles
#echo "snakemake --profile workflow/profiles/local/ $@ $outfiles"
