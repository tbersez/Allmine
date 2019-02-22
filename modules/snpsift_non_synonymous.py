# snpSift SNPs filtration
#
#   This module use snpSift (http://snpeff.sourceforge.net/SnpSift.html)
#   to filter your putative variants. Only non synonymous are kept.
#
#   Input:
#       - sample_variants_annotated.vcf
#
#   Output:
#       - sample_non_synm_variants.vcf
#
#   Parameters:
#       No fancy parameters here...

rule snpSift:
    input:
        vcf = config["VAR"] + "{samples}_varscan_filtered_parsed_annotated.vcf"
    output:
        non_syno = config["VAR"] + "non_synonymous_vars/{samples}_non_syn.vcf"
    message: "Getting non synonimous variants from {input.vcf} \n"
    shell:
        """
        java -jar /snpEff/SnpSift.jar filter \
        -f {input.vcf} \
        --addFilter # TODO: DEFINE THE EXP TO GET NON SYNM VARS \
        > {output.non_syno}
        """

# NB: this module constitute the end point of the pipeline
