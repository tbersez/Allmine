# STAR aligner paired end mode, second pass
#
#   This module runs the second path of the STAR aligner 2 path
#   strategy. The goal is to align reads taking in account splice
#   junction found in the fist pass..
#
#   Inputs:
#       - sample_1_trim.fastq.gz
#       - sample_2_trim.fastq.gz
#       - splicing junction files (.tab)
#
#   Output:
#       - aligned reads
#       - logs for follow up and debuging if needed
#
#   Parameters:
#       No fancy parameters needed, only the threads number is specified.

rule star_pe_SP:
        input:
            # fake input
            flag = config["REF"] + "REindexing_done.txt",
            R1 = config["TRIMMED"] + "{samples}_1_trim.fastq.gz",
            R2 = config["TRIMMED"] + "{samples}_2_trim.fastq.gz",
            genomeDir = config["REF"]
        output:
            bam = config["MAP"] + "{samples}_sorted.bam.gz"
        params:
            prefix = config["MAP"] + "SP/" + "{samples}.",
            tmp = config["MAP"] + "SP/" + "{samples}_SP_STAR_TMP",
            bind = config["BIND"],
            cont = config["CONT"]
        benchmark:
            "benchmarks/star_SP/{samples}.tsv"
        message : "Running STAR second pass with {input.R1} and {input.R2}. \n"
        shell:
            """
            singularity exec -B /{params.bind} {params.cont} \
            STAR \
            --runThreadN 10 \
            --genomeDir {input.genomeDir} \
            --readFilesIn {input.R1} {input.R2} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.prefix} \
            --outStd BAM_SortedByCoordinate \
            --outTmpDir {params.tmp} \
            --outFilterScoreMinOverLread 0.3 \
            --outFilterMatchNminOverLread 0.3 \
            --readFilesCommand zcat | gzip --stdout > {output.bam}
            """
