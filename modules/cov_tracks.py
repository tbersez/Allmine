# Compute depth on bam files
#
#   This module is used to compute depth of reading
#   in regions of interest in your bam files.
#
#   Input:
#       - sample.bam
#       - regions.bed
#
#   Output:
#       - coverage_track.tab
#
#   Parameters:
#       No fancy parameters here...

rule cov_track:
    input:
        bams = expand(config["MAP"] + "{samples}_sorted_parsed.bam", samples = config["samples"])
    params:
        bed = config["REGIONS"]
    output:
        cov = "Coverage_Track.tab",
        bind = config["BIND"],
        cont = config["CONT"]
    message: "Computing depth in regions {params.bed} on bams {input.bam} \n"
    shell:
        """
        singularity exec -B {params.bind} {params.cont} \
        samtools depth \
        -a \
        -b {params.bed} \
        {input.bams} > {output.cov}
        """
