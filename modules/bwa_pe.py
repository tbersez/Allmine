#bwa paired end mode

#tpm path to config and cwd for testing
configfile: "../config_se.yaml"
cwd = os.getcwd() + "/"
#######################################

rule all:
    input:
        expand(cwd + config["MAP"] + "{samples}_trim.fastq.gz", samples = config["samples"])
#######################################

rule run_bwa_paired :
    
