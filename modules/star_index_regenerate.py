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
        denovo_SJ = expand(config["MAP"] + "FP/STAR_SJ/" + "{samples}.SJ.out.tab", samples = config["samples"])
    output:
        # Regenerated genome indexes
        flag = config["REF"] + "REindexing_done.txt"
    params:
        geno_dir = config["REF"],
        bind = config["BIND"],
        cont = config["CONT"]
    benchmark:
            "benchmarks/star_reindex/REindexing.tsv"
    message: "RE-Indexing {input.genome} using de novo SJ from the first pass \n"
    shell:
        """
        singularity exec -B {params.bind} {params.cont} STAR \
        --runMode genomeGenerate \
        --runThreadN 10 \
        --genomeDir  {params.geno_dir}\
        --genomeFastaFiles {input.genome} \
        --sjdbFileChrStartEnd {input.denovo_SJ} \
        --sjdbGTFfile {input.ano}
        touch {output.flag}
        """
