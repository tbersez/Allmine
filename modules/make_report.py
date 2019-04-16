# Report generator
#
#   Starts the .py script for report generation
#
#   Inputs:
#       - sample_non_synm_variants.vcf
#
#   Output:
#       - Non_synonymous_variants_summary.tab
#
#   Parameters:
#       None

rule make_report :
    input:
        non_syno = expand(config["VAR"] + "{samples}/{samples}_varscan.avinput.exonic_variant_function", samples = config["samples"])
    output:
        report = 'Non_synonymous_variants_summary.tab'
    shell:
        """
        ./modules/report_generator.py
        """
