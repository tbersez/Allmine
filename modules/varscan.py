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
#       - P-value <= 0.99 (equivalent to no pvalue filter)
#   Those parameters are quite permisive, thats why a filtering step is needed
#   (they can be modified in the varscan_filtering.py script).

rule varscan:
    input:
        bam = config["MAP"] + "{samples}_sorted_parsed.bam"
    params:
        ref = config["REF"] + config["GENOME"]
    output:
        var = config["VAR"] + "{samples}_varscan.vcf"
    threads: config["THREADS"]
    message: "Looking for SNP in {input.bam} with Varscan \n"
    #from mpileup to varscan to save disk space
    shell:
        """
        singularity exec -B /mnt/nas_eic/gafl01/home/gafl/tbersez ~/Allmine/AllMine \
        samtools mpileup \
        -C 50 \
        -A \
        -f {params.ref} \
        {input.bam} | \
        singularity exec -B /mnt/nas_eic/gafl01/home/gafl/tbersez ~/Allmine/AllMine \
        varscan mpileup2snp \
        --p-value 0.99 \
        --min-coverage 5 \
        --output-vcf 1 \
        --min-reads2 2 \
        --min-var-freq 0.01 \
        --min-avg-qual 15 \
        > {output.var}
        """
