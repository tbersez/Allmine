# STAR aligner single end mode, first path
#
#   This module runs the first pass of the STAR aligner 2 pass
#   strategy. The goal is to discover non annotated splice junctions
#   onto the reference genome used. The sum of all de novo junctions
#   detected will be used for the second path.
#
#   Inputs:
#       - sample_trim.fastq.gz
#       - STAR index and reference genomeDir
#
#   Output:
#       - denovo junctions (tab file)
#       - logs for follow up and debuging if needed
#
#   Parameters:
#       No fancy parameters needed, only the threads number is specified.

include : "./star_se_SP.py"

rule star_se_FP:
    input:
        R1 = config["TRIMMED"] + "{samples}_trim.fastq.gz",
        #fake input to force index building
        ano = config["REF"] + "SAindex",
        genomeDir = config["REF"]
    output:
        denovo_SJ = protected(config["MAP"] + "FP/STAR_SJ/" + "{samples}.SJ.out.tab")
    params:
        prefix = config["MAP"] + "FP/STAR_SJ/" + "{samples}.",
        tmp = directory(config["MAP"] + "FP/" + "{samples}_STAR_TMP")
    message : "Running STAR first pass with {input.R1} to get denovo SJ. \n"
    shell :
        """
        singularity exec ~/Allmine/AllMine STAR \
        --runThreadN  2 \
        --genomeDir {input.genomeDir} \
        --readFilesIn {input.R1} \
        --outTmpDir {params.tmp} \
        --outFileNamePrefix {params.prefix}
        """
