rule download_chain_file_hg19tohg38:
    output:
        'results/refs/ucsc/hg19ToHg38.over.chain.gz'
    shell:
        r"""
            wget -O {output} https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
        """
