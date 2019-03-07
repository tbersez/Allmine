# STAR aligner paired end mode, first path
#
#   This module runs the first pass of the STAR aligner 2 pass
#   strategy. The goal is to discover non annotated splice junctions
#   onto the reference genome used. The sum of all de novo junctions
#   detected will be used for the second path.
#
#   Inputs:
#       - sample_1_trim.fastq.gz
#       - sample_2_trim.fastq.gz
#       - STAR index and reference genomeDir
#
#   Output:
#       - denovo junctions (tab file)
#       - logs for follow up and debuging if needed
#
#   Parameters:
#       No fancy parameters needed, only the threads number is specified.

include : "./star_pe_SP.py"

rule star_pe_FP:
    input:
        R1 = config["TRIMMED"] + "{samples}_1_trim.fastq.gz",
        R2 = config["TRIMMED"] + "{samples}_2_trim.fastq.gz",
        #fake input to force index building
        ano = config["REF"] + "SAindex",
        genomeDir = config["REF"]
    output:
        denovo_SJ = config["MAP"] + "FP/STAR_SJ/" + "{samples}.SJ.out.tab"
    params:
        prefix = config["MAP"] + "FP/STAR_SJ/" + "{samples}.",
        tmp = config["MAP"] + "FP/" + "{samples}_STAR_TMP"
    message : "Running STAR first pass with {input.R1} and {input.R2} to get denovo SJ. \n"
    shell :
        """
        singularity exec -B /mnt/nas_eic/gafl01/home/gafl/tbersez \
         ~/Allmine/AllMine STAR --runThreadN 2 \
        --genomeDir {input.genomeDir} \
        --readFilesIn {input.R1} {input.R2} \
        --outTmpDir {params.tmp} \
        --readFilesCommand zcat \
        --outFilterScoreMinOverLread 0.3 \
        --outFilterMatchNminOverLread 0.3 \
        --outFileNamePrefix {params.prefix}
        """
