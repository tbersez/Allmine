configfile: "./config_pe.yaml"
cwd = os.getcwd() + "/"
# modules loading...
include : cwd + "modules/*"
#############################
rule all:
    input:
        expand(config["TRIMMED"] + "{samples}_R1_trim.fastq.gz", samples = config["samples"]) +
        expand(config["TRIMMED"] + "{samples}_R2_trim.fastq.gz", samples = config["samples"])

rule run_fastp:
    input:
        R1 = lambda wildcards: config["samples"][wildcards.samples]["R1"],
        R2 = lambda wildcards: config["samples"][wildcards.samples]["R2"]
    output:
        R1 =  config["TRIMMED"] + "{samples}_R1_trim.fastq.gz",
        R2 =  config["TRIMMED"] + "{samples}_R2_trim.fastq.gz"
    params:
        title = lambda wildcards: config["samples"][wildcards.samples]["name"]
    threads: config["THREADS"]
    message: "Running fastp on files {input.R1} and {input.R2} \n"
    benchmark:
        "benchmarks/{samples}.fastp.benchmark.txt"
    shell:
        """
        fastp \
        -i {input.R1} \
        -I {input.R2} \
        -o {output.R1} \
        -O {output.R2} \
        -R {params.title} \
        -h {params.title}_QC.html \
        -j {params.title}_QC.json \
        """
