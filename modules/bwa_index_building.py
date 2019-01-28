#bwa paired end mode

#tpm path to config and cwd for testing
configfile: "config_pe.yaml"
cwd = os.getcwd() + "/"
#######################################

rule all:
    input:
        cwd + config["REF"] + config["GENOME"] + ".amb",
        cwd + config["REF"] + config["GENOME"] + ".ann",
        cwd + config["REF"] + config["GENOME"] + ".bwt",
        cwd + config["REF"] + config["GENOME"] + ".pac",
        cwd + config["REF"] + config["GENOME"] + ".sa"

rule bwa_index:
    input:
        genome = cwd + config["REF"] + config["GENOME"]
    output:
        cwd + config["REF"] + config["GENOME"] + ".amb",
        cwd + config["REF"] + config["GENOME"] + ".ann",
        cwd + config["REF"] + config["GENOME"] + ".bwt",
        cwd + config["REF"] + config["GENOME"] + ".pac",
        cwd + config["REF"] + config["GENOME"] + ".sa"
    message: "Building BWA index for reference genome {input.genome}\n"
    threads: config["THREADS"]
    shell:
        """
        /home/aa/anaconda3/bin/bwa index\
        {input.genome}
        """
