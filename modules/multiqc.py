# MultiQC, interactive QC summary
#
#   This module uses multiqc to agregate all fastqc reports
#   in one big reports
#
#   Inputs:
#       - samples_fastqc.html
#
#   Output:
#       - global_qc_report.html
#
#   Parameters:
#       None

rule multiQC:
    input:
        fastqc_flag = expand("QC_post_preproc/{samples}.flag", samples = config["samples"])
    output:
        "Global_QC_summary.html"
    shell:
        """
        singularity exec -B /mnt/nas_eic/gafl01/home/gafl/tbersez ~/Allmine/AllMine \
        multiqc \
        --quiet \
        --outdir ./ \
        --filename Global_QC_summary.html \
        --dirs-depth 1 \
        --fullnames \
        --title Global_QC_summary \
        --comment 'Generated during an AllMine run' \
        ./QC_post_preproc
        """
