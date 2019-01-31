# fastp, single end mode :
#
#   This module run fastp (https://doi.org/10.1093/bioinformatics/bty560)
#   on your single end raw fastq files.
#   fastp perfomrs:
#       - adapters dectection and trimming
#       - low quality bases trimming on both 3' and 5' ends (Qscore < 20)
#       - low complexity regions elimination (ex: polyA tails)
#       - QC report generation in .html format
#
#   Input:
#       - sample.fastq
#
#   Output:
#       - sample_1_trim.fastq
#
#   Parameters:
#       fastp default parameters (see fastp manual for more information).
#       Thoses parameters all well suited for most of sequencing data,
#       However you can modify them directly into the script bellow if you
#       feel the need to (for advanced users).

rule run_fastp_single:
    input:
        R1 = lambda wildcards: config["samples"][wildcards.samples]["R1"]
    output:
        R1 = config["TRIMMED"] + "{samples}_trim.fastq.gz",
        html = cwd + "QC_reports/{samples}_QC.html"
    message: "Running fastp on file {input.R1}\n"
    log: config["LOG"] + "fastp/{samples}.log"
    params:
        title = lambda wildcards: config["samples"][wildcards.samples]["name"]
    threads: config["THREADS"]
    shell:
        """
        /home/aa/anaconda3/bin/fastp \
        -i {input.R1} \
        -o {output.R1} \
        -R {params.title} \
        -h {output.html}

        rm -f fastp.json
        """
