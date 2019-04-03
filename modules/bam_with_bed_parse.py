# Bam file parsing with bed file as blueprint
#
#   This module is used to parse your aligned bam files
#   using the bed file containing the regions of interest
#   as blueprint.
#
#   Input:
#       - sample.bam
#       - regions.bed
#
#   Output:
#       - sample_parsed.bam
#
#   Parameters:
#       No fancy parameters here...

rule parse_bam_with_bed:
    input:
        bam = config["MAP"] + "{samples}_sorted.bam.gz"
    output:
        bam = config["MAP"] + "{samples}_sorted_parsed.bam"
    params:
        bed = config["REGIONS"],
        ref = config["REF"] + config["GENOME"],
        bind = config["BIND"],
        cont = config["CONT"]
    benchmark:
        "benchmarks/bam_parse/{samples}.tsv"
    message: "Parsing {input.bam} using the blueprint {params.bed} \n"
    shell:
        """
        gzip --stdout -d {input.bam} | \
        singularity exec -B {params.bind} {params.cont} samtools view \
        -b \
        -L {params.bed} \
        -O BAM \
        - \
        > {output.bam}
        """
