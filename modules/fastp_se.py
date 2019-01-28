#fastp, single end mode

#tpm path to config and cwd for testing
configfile: "../config_se.yaml"
cwd = os.getcwd() + "/"

rule all:
    input:
        expand(cwd + config["TRIMMED"] + "{samples}_trim.fastq.gz", samples = config["samples"])

rule run_fastp_paired:
    input:
        R1 = lambda wildcards: config["samples"][wildcards.samples]["R1"]
    output:
        R1 = cwd + config["TRIMMED"] + "{samples}_trim.fastq.gz"
    message: "Running fastp on file {input.R1}\n"
    params:
        title = lambda wildcards: config["samples"][wildcards.samples]["name"]
    threads: config["THREADS"]
    shell:
        """
        /home/aa/anaconda3/bin/fastp \
        -i {input.R1} \
        -o {output.R1} \
        -R {params.title} \
        -h cwdreport/{params.title}_QC.html \
        -j cwdreport/{params.title}_QC.json
        """
