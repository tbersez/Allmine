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
#       - Run_report.html
#
#   Parameters:
#       None

rule make_markdown :
    input:
        "Non_synonymous_variants_summary.tab",
        "Coverage_Track.tab"
    output:
        Rmd = "Run_report.html"
    params:
        bind = config["BIND"],
        cont = config["CONT"]
    shell:
        """
        singularity exec -B {params.bind} {params.cont}\
        R -e \
        "rmarkdown::render('./modules/markdown_gen.Rmd', output_file='../{output.Rmd}')"
        """
