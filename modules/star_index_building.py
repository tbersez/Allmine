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

include : "./star_index_regenerate.py"

rule STAR_index:
    input:
        genome = config["REF"] + config["GENOME"],
        ano = config["REF"] + config["ANO"]
    output: # STAR index creates a lot of files!
        config["REF"] + "chrLength.txt",
        config["REF"] + "chrNameLength.txt",
        config["REF"] + "chrName.txt",
        config["REF"] + "chrStart.txt",
        config["REF"] + "exonGeTrInfo.tab",
        config["REF"] + "exonInfo.tab",
        config["REF"] + "geneInfo.tab",
        config["REF"] + "Genome",
        config["REF"] + "genomeParameters.txt",
        config["REF"] + "SA",
        config["REF"] + "SAindex",
        config["REF"] + "sjdbInfo.txt",
        config["REF"] + "sjdbList.fromGTF.out.tab",
        config["REF"] + "sjdbList.out.tab",
        config["REF"] + "transcriptInfo.tab"
    params:
        geno_dir = config["REF"]
    message: "Indexing {input.genome} with {input.ano} for STAR aligner \n"
    shell:
        """
        singularity exec -B /mnt/nas_eic/gafl01/home/gafl/tbersez \
        ~/Allmine/AllMine STAR --runMode genomeGenerate \
        --runThreadN 10 \
        --genomeDir  {params.geno_dir}\
        --genomeFastaFiles {input.genome} \
        --sjdbGTFfile {input.ano}
        mkdir -p mapped/FP/STAR_SJ
        mkdir -p mapped/SP
        """
