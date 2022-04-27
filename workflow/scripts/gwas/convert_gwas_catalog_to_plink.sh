# The main idea is to process GWAS_input_colocalization.GRCh38.txt into input files 
# for qmplot functions.

# tester
input_fn="results/main/2021_Nikhil_eQTL/Data/T1D_GWAS/T1D_34012112_Gaulton/GRCh38/GWAS_input_colocalization.txt"
output_fn="results/main/2021_Nikhil_eQTL/Data/T1D_GWAS/T1D_34012112_Gaulton/GRCh38/GWAS_input_colocalization.qmplot_input.txt"

#input_fn=$1
#output_fn=$2
old_header="CHR	POS	BETA	SE	P	N"
new_header="#CHROM	POS	ID	REF	ALT	A1	TEST	OBS_CT	BETA	SE	T_STAT	P"
echo $new_header > $output_fn
sed '1d' $input_fn | awk 'BEGIN{OFS="	"}; {print $1, $2, "$1:$2", ".", ".", ".", "ADD", $6, $3, $4, ".", $5}' >> $output_fn
