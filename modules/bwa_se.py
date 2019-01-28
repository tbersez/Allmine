#bwa paired end mode

#tpm path to config and cwd for testing
configfile: "config_pe.yaml"
cwd = os.getcwd() + "/"
#######################################

rule all:
    input:
        expand(config["MAP"] + "{samples}_sorted.bam", samples = config["samples"])
#######################################

rule run_bwa_paired :
    input:
        R1 = expand(config["TRIMMED"] + "{samples}_1_trim.fastq.gz", samples = config["samples"])
    output:
        bam = expand(config["MAP"] + "{samples}_sorted.bam", samples = config["samples"])
    params:
        idxbase = config["REF"] + config["GENOME"]
    message: "Mapping reads {input.R1} to {params.idxbase} using BWA.\n"
    shell:
        """
        /home/aa/anaconda3/bin/bwa mem \
        {params.idxbase} \
        {input.R1} \
        | /usr/bin/samtools view -Sb - \
        | /usr/bin/samtools sort -o {output.bam} -
        """
