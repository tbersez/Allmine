# FastQC quality control
#
#   This module run FastQC on your pair end processed reads to
#   ensure their quality. Several features, such as read
#   quality, GC%, dupplication levels ect., are inspected
#   by FastQC. Report are generated in .html format.
#   You can inspect them using a web browser.
#
#   Input:
#       - samples_trim_R1.fastq
#       - samples_trim_R1.fastq
#
#   Output:
#       - samples_report.html
#
#   Parameters:
#       No fancy parameters here ...

rule fastqc_paired:
    input:
        R1 = config["TRIMMED"] + "{samples}_1_trim.fastq.gz",
        R2 = config["TRIMMED"] + "{samples}_2_trim.fastq.gz"
    output:
        dir = "QC_post_preproc/{samples}",
        flag = "QC_post_preproc/{samples}.flag"
    message: "QC on trimmed reads {input.R1} & {input.R2}  with FastQC \n"
    shell:
        """
        mkdir -p {output.dir}
        touch {output.flag}
        singularity exec -B /mnt/nas_eic/gafl01/home/gafl/tbersez ~/Allmine/AllMine \
        fastqc -q \
        --noextract \
        -o {output.dir} \
        {input.R1} {input.R2}
        """
