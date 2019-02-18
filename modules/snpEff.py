# snpEff SNPs annotation
#
#   This module use snpEff (http://snpeff.sourceforge.net/index.html)
#   to annotate your putative variants. Their supposed effect on the
#   coded protein is also computed. For a quick overview of the
#   results, please consult the annotation report (on per sample).
#
#   Input:
#       - sample_varscan_filtered_parsed.vcf
#       - annotation_snpEff_db.tab
#
#   Output:
#       - sample_varscan_filtered_parsed_annotated.vcf
#       - annotation_report.html
#
#   Parameters:
#       Since the purpose of this pipeline is to perform allele mining
#       only UTR and exon variants will be reported. Still, mutation
#       leading to loss of function are also outputed but may not be
#       considered usefull for gene transfert !

rule snpEff:
    input:
        vcf = config["VAR"] + "{samples}_varscan_filtered.vcf"
    output:
        vcf = protected(config["VAR"] + "{samples}_varscan_filtered_parsed_annotated.vcf"),
        rep = protected(config["VAR"] + '{samples}_annotated_report.html')
    message: "Annotating {input.vcf} with snpEff \n"
    params:
        threads = config["THREADS"],
    shell:
        """
        java -jar snpEff.jar \
        eff \
        -i vcf \
        -o vcf \
        -htmlStats {output.rep} \
        -no-downstream \
        -no-intergenic \
        -no-intron \
        -no-upstream \
        -noLog \
        -config ./snpEff.config \
        snpEff_db \
        {input.vcf} > {output.vcf}
        """
