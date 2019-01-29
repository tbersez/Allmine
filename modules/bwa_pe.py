#bwa paired end mode

rule run_bwa_paired :
    input:
        R1 = expand(config["TRIMMED"] + "{samples}_1_trim.fastq.gz", samples = config["samples"]),
        R2 = expand(config["TRIMMED"] + "{samples}_2_trim.fastq.gz", samples = config["samples"]),
        #fake input used to force index building before alignement
        idx = config["REF"] + config["GENOME"] + ".bwt"
    output:
        bam = protected(expand(config["MAP"] + "{samples}_sorted.bam", samples = config["samples"]))
    params:
        idxbase = config["REF"] + config["GENOME"]
    message: "Mapping reads {input.R1} and {input.R2} to {params.idxbase} using BWA.\n"
    #converting to bam, sorting and removing dupplicates in a single command!
    shell:
        """
        /home/aa/anaconda3/bin/bwa mem \
        {params.idxbase} \
        {input.R1} \
        {input.R2} \
        | /usr/bin/samtools view -Sb - \
        | /usr/bin/samtools sort -o - \
        | /usr/bin/samtools rmdup -s - {output.bam}
        """
