# snpSift SNPs filtration
#
#   This module use snpSift (http://snpeff.sourceforge.net/SnpSift.html)
#   to filter your putative variants. Only non synonymous SNPs are kept.
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
        vcf = config["VAR"] + "{samples}/{samples}_varscan_filtered_parsed_annotated.vcf"
    output:
        non_syno = config["VAR"] + "{samples}/{samples}_non_synonymous.vcf"
    message: "Getting non synonymous variants from {input.vcf} \n"
    shell:
        """
        cat {input.vcf} | \
        singularity exec -B /mnt/nas_eic/gafl01/home/gafl/tbersez \
        ~/Allmine/AllMine \
        java -jar /snpEff/SnpSift.jar filter "ANN[*].EFFECT has 'missense_variant'" \
        > {output.non_syno}
        """

# NB: this module constitute the end point of the pipeline
