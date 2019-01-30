# BWA single end mode:
#
#   This use bwa mem (http://bio-bwa.sourceforge.net/) to align
#   your processed single end reads to your genome/transciptome/cds.
#   To save computation time, pipes are use to convert the bwa output
#   to the bam format, sort them and remove PCR dupplicates using
#   Samtools (http://samtools.sourceforge.net/).
#
#   Input:
#       - sample_trim.fastq
#
#   Output:
#       - sample_sorted.bam
#
#   Parameters:
#       bwa mem default parameters (see bwa manual for details).
#       Advanced users can modify them in the script bellow if needed.
#       Note : the bwa mem algorithm was choosed among others bwa algorithms
#       because recomanded for s generally recommended for high-quality
#       queries as it is faster and more accurate.

rule run_bwa_paired :
    input:
        R1 = expand(config["TRIMMED"] + "{samples}_1_trim.fastq.gz", samples = config["samples"]),
        #fake input used to force index building before alignement
        idx = config["REF"] + config["GENOME"] + ".bwt"
    output:
        bam = protected(expand(config["MAP"] + "{samples}_sorted.bam", samples = config["samples"]))
    params:
        idxbase = config["REF"] + config["GENOME"]
    message: "Mapping reads {input.R1} to {params.idxbase} using BWA.\n"
    #converting to bam, sorting and removing dupplicates in a single command!
    shell:
        """
        /home/aa/anaconda3/bin/bwa mem \
        {params.idxbase} \
        {input.R1} \
        | /usr/bin/samtools view -Sb - \
        | /usr/bin/samtools sort -o - \
        | /usr/bin/samtools rmdup -s - {output.bam}
        """
