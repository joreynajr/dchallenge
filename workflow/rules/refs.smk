# rule to download the chain file converting between hg19 and hg38
rule download_chain_file_hg19tohg38: #(status: running)
    output:
        'results/refs/ucsc/hg19ToHg38.over.chain.gz'
    shell:
        r"""
            wget -O {output} https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
        """

# rule to download the chain file converting between hg38 and hg19
# currently the link is wrong.
#rule download_chain_file_hg38tohg19: #(status: developing)
#    output:
#        'results/refs/ucsc/hg38ToHg19.over.chain.gz'
#    shell:
#        r"""
#            wget -O {output} https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
#        """
