# dchallenge
As the Ay Lab we are competing in the DChallenge where we are looking SNPs colocalizing for T1D and expression and integrating this with HiChIP interactions. 

We are futher documenting the pipeline using protocols.io for which you can using the following url:
https://www.protocols.io/private/B58E0187427B11ECAA1D0A58A9FEAC02

Lastly, we are making a submission write up. Only those with permission can visit this file at this moment:
https://docs.google.com/document/d/1ecb7tHfzQIndU_BC9SuwtD-95neMvxdVKeP-zvjhbBc/edit?usp=sharing

# Pipeline
## Analyzing loops with colocalized SNP-gene pairs 
1) Downloaded T1D GWAS data for:
    - (Chiou et al 2021) I commonly reference it as Gaulton data since he is the PI
        link: https://www.nature.com/articles/s41586-021-03552-w
    - (tbd) I am working on finding the Stephen data
2) Downloaded eQTL data from the (Mu et al., 2021) paper
    link: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02334-x
3) Ran the colocalization analysis pipeline prepared by Sourya. Uses the coloc package from R. 
    link: https://cran.r-project.org/web/packages/coloc/index.html
4) Got FitHiChIP results from Sourya
5) Downloaded the GENCODE v19 gene annotation file (GFF)
6) Investigating the connections between SNPs, genes and loops (denoting these SGLoops)
7) Visualized these results in the WashU Browser
8) Further investigated genes in SGLoops to find high confidence gene targets
