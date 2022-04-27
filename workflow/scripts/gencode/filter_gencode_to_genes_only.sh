input="results/refs/gencode/v30/gencode.v30.annotation.grch37.bed"
output="results/refs/gencode/v30/gencode.v30.annotation.grch37.genes_only.bed"
read -r cmd << EOM
awk 'BEGIN{OFS="\t"}; {if (\$5 == "gene") print \$0};' $input > $output
EOM
echo "Running: $cmd"
