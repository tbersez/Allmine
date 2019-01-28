#fastp, paired end mode

#tpm path to config and cwd for testing
configfile: "../config_pe.yaml"
cwd = os.getcwd() + "/"
######################################

rule all:
    input:
        expand(cwd + config["TRIMMED"] + "{samples}_1_trim.fastq.gz", samples = config["samples"]) +
        expand(cwd + config["TRIMMED"] + "{samples}_2_trim.fastq.gz", samples = config["samples"])
######################################

rule run_fastp_paired:
    input:
        R1 = lambda wildcards: config["samples"][wildcards.samples]["R1"],
        R2 = lambda wildcards: config["samples"][wildcards.samples]["R2"]
    output:
        R1 = cwd + config["TRIMMED"] + "{samples}_1_trim.fastq.gz",
        R2 = cwd + config["TRIMMED"] + "{samples}_2_trim.fastq.gz",
        html = cwd + "QC_reports/{samples}_QC.html"
    message: "Running fastp on files {input.R1} and {input.R2} \n"
    params:
        title = lambda wildcards: config["samples"][wildcards.samples]["name"]

    threads: config["THREADS"]
    shell:
        """
        /home/aa/anaconda3/bin/fastp \
        -i {input.R1} \
        -I {input.R1} \
        -o {output.R1} \
        -O {output.R2} \
        -R {params.title} \
        -h {output.html}\

        rm -f fastp.json
        """
