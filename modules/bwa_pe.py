# BWA paired end mode:
#
#   This module use bwa mem (http://bio-bwa.sourceforge.net/) to align
#   your processed paired end reads to your genome/transciptome/cds.
#   To save computation time, pipes are use to convert the bwa output
#   to the bam format, sort them and remove PCR dupplicates using
#   Samtools (http://samtools.sourceforge.net/).
#
#   Input:
#       - sample_1_trim.fastq
#       - sample_2_trim.fastq
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

rule bwa_paired :
    input:
        R1 = config["TRIMMED"] + "{samples}_1_trim.fastq.gz",
        R2 = config["TRIMMED"] + "{samples}_2_trim.fastq.gz",
        #fake input used to force index building before alignement
        idx = config["REF"] + config["GENOME"] + ".bwt"
    output:
        bam = temp(config["MAP"] + "{samples}_sorted.bam"),
        bam_gz = config["MAP"] + "{samples}_sorted.bam.gz"
    params:
        idxbase = config["REF"] + config["GENOME"],
        bind = config["BIND"],
        cont = config["CONT"]
    benchmark:
        "benchmarks/bwa/{samples}.tsv"
    message: "Mapping reads {input.R1} and {input.R2} to {params.idxbase} using BWA.\n"
    #converting to bam, sorting and removing dupplicates in a single command!
    shell:
        """
        singularity exec -B {params.bind} {params.cont} bwa mem \
        -t 10 \
        -w 100 \
        -d 100 \
        -r 1.5 \
        -B 4 \
        -O 6 \
        -E 1 \
        -L 5 \
        {params.idxbase} \
        {input.R1} \
        {input.R2} \
        | /usr/bin/samtools view -Sb -@ 10 - \
        | /usr/bin/samtools sort -@ 10 -o - \
        | /usr/bin/samtools rmdup -s - {output.bam}
        gzip --keep {output.bam}
        """
