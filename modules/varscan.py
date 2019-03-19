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
#       - Minimum read depth at a position
#       - Minimum base quality at a position
#       - Minimum variant allele frequency threshold
#       - P-value <= 0.99 (equivalent to no pvalue filter)
#   Thoses parameters may have to be edited to tune the calling sensitivity !

rule varscan:
    input:
        bam = config["MAP"] + "{samples}_sorted_parsed.bam"
    params:
        ref = config["REF"] + config["GENOME"],
        bind = config["BIND"],
        cont = config["CONT"]
    output:
        var = config["VAR"] + "{samples}/{samples}_varscan.vcf",
    message: "Looking for SNP in {input.bam} with Varscan \n"
    benchmark:
        "benchmarks/varscan/{samples}.tsv"
    #from mpileup to varscan to save disk space
    shell:
        """
        singularity exec -B {params.bind} {params.cont} \
        samtools mpileup \
        -C 50 \
        -A \
        -f {params.ref} \
        {input.bam} | \
        singularity exec -B {params.bind} {params.cont} \
        varscan mpileup2snp \
        --p-value 0.99 \
        --min-coverage 8 \
        --output-vcf 1 \
        --min-strands2 1 \
        --min-var-freq 0.05 \
        --min-avg-qual 20 \
        > {output.var}
        """
