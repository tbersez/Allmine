# VarScan variant calling
#
#   This module run VarScan (http://varscan.sourceforge.net/)
#   on your aligned reads (e.i. bam files) to detect SNP.
#   VarScan takes a pileup file as input, which is piped in
#   using Samtools (http://samtools.sourceforge.net/).
#
#   Input:
#       - sample_sorted.bam
#
#   Output:
#       - sample_varscan.tab
#
#   Parameters:
#     A putative SNP is called if:
#       - Minimum read depth at a position >= 8
#       - Minimum supporting reads at a position >= 2
#       - Minimum base quality at a position >= 15 (Qscore)
#       - Minimum variant allele frequency threshold >= 0.01
#       - P-value <= 0.05
#   Those parameters are quite permisive, thats why a filtering step is needed
#   (they can be modified in the varscan_filtering.py script).

rule run_varscan:
    input:
        bam = expand(config["MAP"] + "{samples}_sorted.bam", samples = config["samples"])
    params:
        ref = config["REF"] + config["GENOME"]
    output:
        var = protected(expand(config["VAR"] + "{samples}_varscan.tab", samples = config["samples"]))
    threads: config["THREADS"]
    message: "Looking for SNP & indels in {input.bam} with Varscan \n"
    #from mpileup to varscan to save disk space
    shell:
        """
        /usr/bin/samtools mpileup -A \
        -f {params.ref} \
        {input.bam} \
        | /home/aa/anaconda3/bin/varscan pileup2snp \
        --p-value 0.05 - \
        > {output.var}
        """