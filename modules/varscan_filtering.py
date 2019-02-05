# VarScan filtering variants
#
#   This module run VarScan filter (http://varscan.sourceforge.net/)
#   on your putative SNPs. Filtering reduce false dicovery risks.
#
#   Input:
#       - sample_varscan.tab
#
#   Output:
#       - sample_varscan_filtered.tab
#
#   Parameters:
#     Thoses filters are applied:
#       - Minimum read depth at a position >= 10
#       - Minimum supporting reads at a position >= 2
#       - Minimum base quality at a position >= 20 (Qscore)
#       - Minimum variant allele frequency threshold >= 0.20
#       - P-value <= 0.05
#   Those parameters may be modified in the script below to tune
#   stringency of the filtering.

rule filter_varscan:
    input:
        var = config["VAR"] + "{samples}_varscan.vcf"
    output:
        var = protected(config["VAR"] + "{samples}_varscan_filtered.vcf")
    message: "Filtering Varscan variants from {input.var} \n"
    log: config["LOG"] + "varscan_filtering/{samples}.log"
    threads: config["THREADS"]
    # Appling filter function from Varscan, parameters may be changed to fit your needs
    shell:
        """
        /home/aa/anaconda3/bin/varscan filter \
        {input.var}  \
        --min-coverage 0 \
        --min-reads2 10 \
        --min-strands2 1 \
        --min-avg-qual 20 \
        --min-var-freq 0.20 \
        --p-value 0.05 \
        > {output.var}
        """
