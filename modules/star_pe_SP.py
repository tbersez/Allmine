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
            genomeDir = config["REF"],
            denovo_SJ = expand(config["MAP"] + "STAR_SJ/" + "{samples}.SJ.out.tab", samples = config["samples"])
        output:
            bam = config["MAP"] + "{samples}_sorted.bam"
        params:
            prefix = config["MAP"] + "{samples}.",
            threads = config["THREADS"],
            tmp = config["MAP"] + "{samples}_sp_STAR_TMP"
        log: config["LOG"] + "STAR_SP/{samples}.log"
        message : "Running STAR second pass with {input.R1} and {input.R2}. \n"
        shell:
            """
            STAR \
            --runThreadN {params.threads} \
            --genomeDir {input.genomeDir} \
            --readFilesIn {input.R1} {input.R2} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.prefix} \
            --outStd  BAM_SortedByCoordinate \
            --outTmpDir {params.tmp} \
            --readFilesCommand zcat > {output.bam}
            """
