#########################################################################################
# Getting all eQTL datasets
#########################################################################################

## setting the gwas source 
##gwas_source="results/main/2021_Nikhil_eQTL/Data/T1D_GWAS/T1D_34012112_Gaulton/GWAS_input_colocalization_pval_lt_5eMinus8.txt"
#gwas_source="T1D_34012112_Gaulton"
#
#skip_outputs=("results/main/coloc/Results/Colocalization_SMKN/T1D_34012112_Gaulton/GENCORD/LCL/" \
#    "results/main/coloc/Results/Colocalization_SMKN/T1D_34012112_Gaulton/Schmiedel_2018/NK-cell_naive/")
#
## running through all the eqtl files
#outfiles+=""
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
#    new_fn="results/main/coloc/Results/Colocalization_SMKN/${gwas_source}/${eqtl_source}/${ge_source}/"
#    if [[ ! " ${skip_outputs[*]} " =~ " ${new_fn} " ]];
#    then
#        outfiles+=" $new_fn"
#        echo "$fn ----------> $new_fn"
#    fi
#done
#

fns=""
while IFN= read -r line;
do
    run_info=($line)
    fn="results/main/coloc/Results/eQTL_Catalogue/${run_info[0]}/${run_info[1]}/${run_info[2]}/"
    fns+="$fn "

done < <(sed '1d' config/coloc_samplesheets/t1d.samplesheet.tsv | grep -v "^#" | head -n 100)
echo $fns

#outfiles=$(echo $outfiles | cut -d " " -f 4)
#echo
#echo "####################################################################################"
#echo "outfiles: $outfiles"
#echo "####################################################################################"

# TEST Example 
#snakemake --profile workflow/profiles/pbs-torque $@ results/main/coloc/Results/eQTL_Catalogue/T1D_34012112_Gaulton/Schmiedel_2018/test
#snakemake --profile workflow/profiles/local $@ results/main/coloc/Results/eQTL_Catalogue/T1D_34012112_Gaulton/Schmiedel_2018/test

# AD HOC re-run of certain cells
#new_fn="results/main/coloc/Results/Colocalization_SMKN/T1D_34012112_Gaulton/Schmiedel_2018/NK-cell_naive/"
#snakemake --profile workflow/profiles/pbs-torque/ $@ $new_fn
#new_fn="results/main/coloc/Results/eQTL_Catalogue/T1D_34012112_Gaulton/Schmiedel_2018/CD8_T-cell_anti-CD3-CD28"
#snakemake --profile workflow/profiles/pbs-torque/ $@ $new_fn
#fns="results/main/coloc/Results/eQTL_Catalogue/T1D_34012112_Gaulton/Schmiedel_2018/CD8_T-cell_anti-CD3-CD28/"

# FULL Run
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
#snakemake --profile workflow/profiles/local/ $@ $fns
