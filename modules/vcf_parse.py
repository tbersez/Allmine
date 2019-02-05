# Vcftools parsing vcf with bed as blueprint
#
#   This module run vcftools on your putative SNPs in order
#   to extract variants in regions of interest.
#
#   Input:
#       - sample_varscan_filtered.vcf
#       - regions_of_interest.vcf
#
#   Output:
#       - sample_varscan_filtered_parsed.vcf
#
#   Parameters:
#       No fancy parameters here...

rule vcf_bed_parse:
    input:
        var = config["VAR"] + "{samples}_varscan_filtered.vcf"
    output:
        var_bed = config["VAR"] + "{samples}_varscan_filtered_parsed.vcf"
    params:
        bed = config["REGIONS"]
    message: "Parsing {input.var} using {params.bed} with vcftools. \n"
    threads: config["THREADS"]
    shell:
        """
        vcftools \
        --vcf {input.var} \
        --bed {params.bed} \
        --out {output.var_bed} \
        --recode \
        --keep-INFO-all
        """
