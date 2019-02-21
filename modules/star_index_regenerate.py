# STAR index building:
#
#   This module build the STAR (https://github.com/alexdobin/STAR)
#   alignment index from your reference genome and the gff/gtf annotation
#   file associated. This step is computationnaly costly!
#   It must not be runned on a desktop.
#
#   Input:
#       - genome
#       - genome annotation
#
#   Output:
#       - lot of files ! (see bellow)
#       All thoses files are created in the genome directory.
#       They MUST NOT be moved, modified or renamed!
#
#   Parameters:
#       Base command for STAR indexing does not ask for fancy parameters.

rule STAR_REindex:
    input:
        genome = config["REF"] + config["GENOME"],
        ano = config["REF"] + config["ANO"],
        denovo_SJ = expand(config["MAP"] + "STAR_SJ/" + "{samples}.SJ.out.tab", samples = config["samples"])
    output:
        # Regenerated genome indexes
        flag = config["REF"] + "REindexing_done.txt"
    params:
        threads = config["THREADS"],
        geno_dir = config["REF"]
    message: "RE-Indexing {input.genome} using de novo SJ from the first pass \n"
    shell:
        """
        touch {output.flag}
        STAR --runMode genomeGenerate \
        --runThreadN {params.threads} \
        --genomeDir  {params.geno_dir}\
        --genomeFastaFiles {input.genome} \
        --sjdbFileChrStartEnd {input.denovo_SJ}
        """
