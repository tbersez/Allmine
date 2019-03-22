# FastQC quality control
#
#   This module run FastQC on your single end processed reads to
#   ensure their quality. Several features, such as read
#   quality, GC%, dupplication levels ect., are inspected
#   by FastQC. Report are generated in .html format.
#   You can inspect them using a web browser.
#
#   Input:
#       - samples_trim.fastq
#
#   Output:
#       - samples_report.html
#
#   Parameters:
#       No fancy parameters here ...

rule fastqc_single:
    input:
        R1 = config["TRIMMED"] + "{samples}_trim.fastq.gz"
    output:
        dir = 'QC_post_preproc/{samples}',
        flag = "QC_post_preproc/{samples}.flag"
    params:
        bind = config["BIND"],
        cont = config["CONT"]
    benchmark:
        "benchmarks/fastqc/{samples}.tsv"
    message: "QC on trimmed reads {input.R1} with FastQC \n"
    shell:
        """
        mkdir -p {output.dir}
        touch {output.flag}
        singularity exec -B {params.bind} {params.cont} \
        fastqc -q \
        --noextract \
        -o {output.dir} \
        {input.R1}
        """
