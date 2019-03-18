# fastp, paired end mode :
#
#   This module run fastp (https://doi.org/10.1093/bioinformatics/bty560)
#   on your paired end raw fastq files.
#   fastp perfomrs:
#       - adapters dectection and trimming
#       - low quality bases trimming on both 3' and 5' ends (Qscore < 20)
#       - low complexity regions elimination (ex: polyA tails)
#       - QC report generation in .html format
#
#   Input:
#       - sample_1.fastq
#       - sample_2.fastq
#
#   Output:
#       - sample_1_trim.fastq
#       - sample_2_trim.fastq
#
#   Parameters:
#       fastp default parameters (see fastp manual for more information).
#       Thoses parameters all well suited for most of sequencing data,
#       However you can modify them directly into the script bellow if you
#       feel the need to (for advanced users).

rule fastp_paired:
    input:
        R1 = lambda wildcards: config["samples"][wildcards.samples]["R1"],
        R2 = lambda wildcards: config["samples"][wildcards.samples]["R2"]
    output:
        R1 = config["TRIMMED"] + "{samples}_1_trim.fastq.gz",
        R2 = config["TRIMMED"] + "{samples}_2_trim.fastq.gz",
    message: "Running fastp on files {input.R1} and {input.R2} \n"
    params:
        title = lambda wildcards: config["samples"][wildcards.samples]["name"],
        html = config["TRIMMED"] + "{samples}_out.html",
        json = config["TRIMMED"] + "{samples}_out.json",
        bind = config["BIND"],
        cont = config["CONT"]
    shell:
        """
        singularity exec -B {params.bind} {params.cont} fastp \
        -i {input.R1} \
        -I {input.R2} \
        -o {output.R1} \
        -O {output.R2} \
        -R {params.title} \
        -h {params.html} \
        -j {params.json} \
        -w 10

        rm -f {params.html} {params.json}
        """
