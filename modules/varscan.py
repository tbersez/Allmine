#varscan variant calling

#tpm path to config and cwd for testing
configfile: "config_pe.yaml"
cwd = os.getcwd() + "/"
#######################################

rule all:
    input:
        expand(config["VAR"] + "{samples}_varscan.tab", samples = config["samples"])
#######################################

rule run_varscan:
    input:
        bam = expand(config["MAP"] + "{samples}_sorted.bam", samples = config["samples"])
    params:
        ref = config["REF"] + config["GENOME"]
    output:
        var = protected(expand(config["VAR"] + "{samples}_varscan.tab", samples = config["samples"]))
    threads: config["THREADS"]
    message: "Looking for SNP & indels in {input.bam} with Varscan \n"
    #from mpileup to varscan to save disk space
    shell:
        """
        /usr/bin/samtools mpileup \
        -f {params.ref} \
        {input.bam} \
        | /home/aa/anaconda3/bin/varscan pileup2snp
        --p-value 0.05 - \
        > {output.var}
        """
