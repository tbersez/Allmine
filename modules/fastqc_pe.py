# FastQC quality control
#
#   This module run FastQC on your pair end processed reads to
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

rule fastqc:
    input:
        R1 = config["TRIMMED"] + "{samples}_1_trim.fastq.gz",
        R2 = config["TRIMMED"] + "{samples}_2_trim.fastq.gz"
    output:
        dir = directory('QC_post_preproc/{samples}')
    message: "QC on trimmed reads {input.R1} & {input.R2}  with FastQC \n"
    threads: config["THREADS"]
    shell:
        """
        mkdir {output.dir}
        fastqc -o {output.dir} {input.R1} {input.R2}
        """
