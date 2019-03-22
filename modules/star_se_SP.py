# STAR aligner single end mode, second pass
#
#   This module runs the second pass of the STAR aligner 2 path
#   strategy. The goal is to align reads taking in account splice
#   junction found in the fist pass..
#
#   Inputs:
#       - sample_trim.fastq.gz
#       - splicing junction files (.tab)
#
#   Output:
#       - aligned reads
#       - logs for follow up and debuging if needed
#
#   Parameters:
#       No fancy parameters needed, only the threads number is specified.

rule star_se_SP:
        input:
            # fake input
            flag = config["REF"] + "REindexing_done.txt",
            R1 = config["TRIMMED"] + "{samples}_trim.fastq.gz",
            genomeDir = config["REF"]
        output:
            bam = config["MAP"] + "{samples}_sorted.bam"
        params:
            prefix = config["MAP"] + "{samples}.",
            tmp = config["MAP"] + "SP/" + "{samples}_sp_STAR_TMP",
            bind = config["BIND"],
            cont = config["CONT"]
        benchmark:
            "benchmarks/star_SP/{samples}.tsv"
        message : "Running STAR second pass with {input.R1}. \n"
        shell:
            """
            singularity exec -B {params.bind} {params.cont} \
            STAR \
            --runThreadN 10 \
            --genomeDir {input.genomeDir} \
            --readFilesIn {input.R1} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.prefix} \
            --outStd  BAM_SortedByCoordinate \
            --outTmpDir {params.tmp} \
            --readFilesCommand zcat > {output.bam}
            """
