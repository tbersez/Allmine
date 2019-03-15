# BWA index building:
#
#   This module build the bwa (http://bio-bwa.sourceforge.net/)
#   alignment index from your reference genome/transciptome/cds.
#   For large genomes (e.i mamals or plants) the indexation is
#   computationnaly costly, and so, may no be runned on low
#   resources.
#
#   Input:
#       - genome/transciptome/cds.fa
#
#   Output:
#       - genome/transciptome/cds.fa.amb
#       - genome/transciptome/cds.fa.ann
#       - genome/transciptome/cds.fa.bwt
#       - genome/transciptome/cds.fa.pac
#       - genome/transciptome/cds.fa.sa
#       All thoses files are created in the genome directory.
#       They MUST NOT be moved, modified or renamed!
#
#   Parameters:
#       bwa index defaults parameters (see bwa manual for more information).
#       Default building parameters are appropriate in 99% of cases. Advanced
#       users can modify them in the script bellow if needed.

rule bwa_index:
    input:
        genome = config["REF"] + config["GENOME"]
    output:
        protected(config["REF"] + config["GENOME"] + ".amb"),
        protected(config["REF"] + config["GENOME"] + ".ann"),
        protected(config["REF"] + config["GENOME"] + ".bwt"),
        protected(config["REF"] + config["GENOME"] + ".pac"),
        protected(config["REF"] + config["GENOME"] + ".sa")
    message: "Building BWA index for reference genome {input.genome}\n"
    shell:
        """
        singularity exec -B /mnt/nas_eic/gafl01/home/gafl/tbersez ~/Allmine/AllMine bwa index \
        {input.genome}
        """
