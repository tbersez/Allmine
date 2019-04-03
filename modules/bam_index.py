# Samtools index bam parsed files
#
#   This module build index form bam files used by
#   WhatsHap for phasing
#
#   Input:
#       - parsed.bam
#
#   Output:
#       - parsed.bam.fai
#
#   Parameters:
#       None
rule bam_index:
    input:
        bam = config["MAP"] + "{samples}_sorted_parsed.bam"
    output:
        index = config["MAP"] + "{samples}_sorted_parsed.bam.bai"
    params:
        bind = config["BIND"],
        cont = config["CONT"]
    benchmark:
        "benchmarks/samtools/{samples}_index.tsv"
    message: "indexing parsed bam file {input.bam} using Samtools."
    shell:
        """
        singularity exec -B {params.bind} {params.cont} \
        samtools index -b {input.bam} \
        {output.index}
        """
