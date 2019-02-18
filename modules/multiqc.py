#multiQC

#tpm path to config and cwd for testing
configfile: "../config_se.yaml"
cwd = os.getcwd() + "/"
#######################################

rule all:
    input:
        cwd + "QC_reports/Global_QC_summary.html"
#######################################

rule run_multiQC:
    output:
        cwd + "QC_reports/Global_QC_summary.html"
    shell:
        """
        multiqc \
        --quiet \
        --outdir {cwd}QC_reports/ \
        --filename Global_QC_summary.html \
        {cwd}QC_reports
        """
# TODO: fix the issue with fastp/multiQC versions
# for the moment, run the pipe with the -k flag
