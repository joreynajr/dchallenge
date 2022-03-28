# File to describe different modules which are also
# listed in order.

# Module to do some basic rule testing
test_rules.smk

# Basic wrapper rules which can be inherited
wrappers.smk

# Reference file downloading and processing
refs.smk

# Colocalization modules to process GWAS and eQTL data.
# Currently GWAS is downloaded manually.
gwas.smk
coloc.smk

# Finemapping module 
finemapping.smk

# pieQTL module 
pieqtls.smk

# Locating SNP-Gene Loops
sgls.smk

# Utilize SCAVENGE with Chiou et al., 2021 data
scavenge.smk
