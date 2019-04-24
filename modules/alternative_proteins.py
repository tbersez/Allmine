# Alternative proteins generator
#
#   This module build alternative proteins based on
#   non synonymous SNPs annotated by Annovar.
#
#   Input:
#       - {samples}/{samples}_varscan.avinput.exonic_variant_function
#
#   Output:
#       - mutated_proteins.fasta
#
#   Parameters:
#       None...

rule alternative_proteins:
    input:
        exonic_vars = config["VAR"] + "{samples}/{samples}_varscan.avinput.exonic_variant_function"
    output:
        prots = config["VAR"] + "{samples}/{samples}_mutated_proteins.fasta"
    params:
        genome = config["REF"] + config["GENOME"],
        bind = config["BIND"],
        cont = config["CONT"]
    message: "Writing mutated proteins from {input.exonic_vars} variants \n"
    shell:
        """
        singularity exec -B {params.bind} {params.cont} \
        /opt/annovar/coding_change.pl \
        --includesnp \
        --onlyAltering \
        {input.exonic_vars} \
        avdb/AV_refGene.txt \
        avdb/AV_refGeneMrna.fa \
        > {output.prots}
        """
