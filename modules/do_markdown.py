# Markdown builder
#
#   Starts the .R script for Markdown generation
#
#   Inputs:
#       - Coverage_Track.tab
#       - regions.bed
#       - Non_synonymous_variants_summary.tab
#
#   Output:
#       -
#
#   Parameters:
#       None

rule make_report :
    input:
        "Non_synonymous_variants_summary.tab",
        "Coverage_Track.tab"
    output:
        Rmd = "Run_report.html"
    shell:
        """
        R -e \
        "rmarkdown::render('./modules/markdown_gen.Rmd', output_file='../{output.Rmd}')"
        """
