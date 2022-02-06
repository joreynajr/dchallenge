#########################################################################################
# Getting all eQTL datasets
#########################################################################################

# setting the gwas source 
#gwas_source="results/main/2021_Nikhil_eQTL/Data/T1D_GWAS/T1D_34012112_Gaulton/GWAS_input_colocalization_pval_lt_5eMinus8.txt"
gwas_source="T1D_34012112_Gaulton"
outfiles+=""

########################################################################################
# eQTL Catalog files
########################################################################################
# running through all the eqtl files
#for fn in $(ls results/main/eqtl/*/ge/*_ge_*.all.tsv.gz);
#do 
#    # extracting the eqtl_source and ge_source
#    keynames=$(echo $fn | xargs basename | sed -e "s/_ge_/     /" -e "s/\.all\.tsv\.gz//")
#    eqtl_source=$(echo $keynames | cut -d " " -f 1)
#    ge_source=$(echo $keynames | cut -d " " -f 2 | cut -d . -f 1)
#
#    #echo "keynames: $keynames"
#    #echo "eqtl_source: $eqtl_source"
#    #echo "ge_source: $ge_source"
#    #echo
#
#    # getting the output files
#    new_fn="results/main/coloc/Results/Colocalization_SMKN_FDR/${gwas_source}/${eqtl_source}/${ge_source}/"
#    if [ $new_fn != "results/main/coloc/Results/Colocalization_SMKN/T1D_34012112_Gaulton/GENCORD/LCL/" ];
#    then
#        outfiles+=" $new_fn"
#        echo "$fn ----------> $new_fn"
#    fi
#done

#########################################################################################
# ImmuNexUT Catalog files
#########################################################################################
# running through all the eqtl files

eqtl_source="ImmuNexUT"
outfiles+=""
for fn in $(ls results/main/eqtl/ImmuNexUT/raw/E-GEAD-420/*_nominal.txt);
do 
    ge_source=$(basename $fn | sed 's/_nominal\.txt//')

    # getting the output files
    new_fn="results/main/coloc/Results/Colocalization_SMKN_FDR_ImmuNexUT/${gwas_source}/${eqtl_source}/${ge_source}/"
    outfiles+=" $new_fn"
    echo "$fn ----------> $new_fn"
    break
done
#
#
##outfiles=$(echo $outfiles | cut -d " " -f 4)
#echo
#echo "####################################################################################"
##outfiles="results/main/coloc/Results/Colocalization_SMKN_FDR/T1D_34012112_Gaulton/BLUEPRINT/monocyte/"
#echo "outfiles: $outfiles"
#echo "####################################################################################"


##snakemake --profile workflow/profiles/pbs-torque/ -n $outfiles
snakemake --profile workflow/profiles/pbs-torque/ $@ $outfiles
#snakemake --profile workflow/profiles/local/ $@ $outfiles
#echo "snakemake --profile workflow/profiles/local/ $@ $outfiles"
