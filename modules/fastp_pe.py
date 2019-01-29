#fastp, paired end mode

rule run_fastp_paired:
    input:
        R1 = lambda wildcards: config["samples"][wildcards.samples]["R1"],
        R2 = lambda wildcards: config["samples"][wildcards.samples]["R2"]
    output:
        R1 = config["TRIMMED"] + "{samples}_1_trim.fastq.gz",
        R2 = config["TRIMMED"] + "{samples}_2_trim.fastq.gz",
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
