# Annovar SNPs annotation
#
#   This module annotate SNPs from Varscan
#   on a gene model basis
#
#   Input:
#       - sample_varscan.avinput
#
#   Output:
#       - sample_varscan_anno.avoutput
#
#   Parameters:
#       None...

rule annovar :
    input:
        inav = config["VAR"] + "{samples}/{samples}_varscan.avinput"
    output:
        all_vars = config["VAR"] + "{samples}/{samples}_varscan.avinput.variant_function",
        exonic_vars = config["VAR"] + "{samples}/{samples}_varscan.avinput.exonic_variant_function"
    params:
        bind = config["BIND"],
        cont = config["CONT"]
    message: "Annotating {input.inav} with Annovar \n"
    shell:
        """
        singularity exec -B {params.bind} {params.cont} \
        /opt/annovar/annotate_variation.pl \
        --geneanno \
        --dbtype refGene \
        --buildver AV \
        {input.inav} \
        avdb/ \
        """
