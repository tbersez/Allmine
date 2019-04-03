# WhatsHap read-bases phasing:
#
#   This runs WhatsHap on your variant files to try to match
#   phasing based on read information. In brief, reads sharing variants
#   are assigned to the same phase. Haplotypes are reconstructed by iteration.
#
#   Input:
#       - variants.vcf
#       - reference.fna
#       - reference.fna.fai
#
#   Output:
#       - phased_variants.vcf
#
#   Parameters:
#       -multi-threading and max coverage used (X = 20x).
#       Warning !: Computation time goes beyond exponantial
#                  for X greater than 20 !!!

rule whatshap:
    input:
        genome = config["REF"] + config["GENOME"],
        index_genome = config["REF"] + config["GENOME"] + ".fai",
        index_bam = config["MAP"] + "{samples}_sorted_parsed.bam.bai",
        bam = config["MAP"] + "{samples}_sorted_parsed.bam",
        var = config["VAR"] + "{samples}/{samples}_varscan.vcf"
    output:
        phased = config["VAR"] + "{samples}/{samples}_varscan_phased.vcf"
    params:
        bind = config["BIND"],
        cont = config["CONT"]
    benchmark:
        "benchmarks/whatshap/{samples}_whatshap.tsv"
    message: "Phasing variant from {input.var}\n using WhatsHap."
    shell:
        """
        singularity exec -B {params.bind} {params.cont} \
        whatshap phase \
        -H 20 \
        --reference {input.genome} \
        -o {output.phased} \
        {input.var} \
        {input.bam}
        """
