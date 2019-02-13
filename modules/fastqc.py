# FastQC quality control
#
#   This module run FastQC on your processed reads to
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
        dir = config["TRIMMED"]
    output: "QC"
    message: "QC on trimmed reads with FastQC \n"
    threads: config["THREADS"]
    shell:
        """
        fastqc -o QC/ {input.dir}
        """
