input_fn="results/main/2021_Nikhil_eQTL/Data/T1D_GWAS/T1D_34012112_Gaulton/GRCh38/GWAS_input_colocalization.qmplot_input.txt"
title="T1D_25751624"
output_prefix="results/gwas/summary/T1D_34012112_Gaulton.all_gwas"
qmplot -I ${input_fn} -T ${title} --dpi 600 -O ${output_prefix}
