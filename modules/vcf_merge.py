# Gather non synonymous SNPs in one VCF file
#
#   This module use awk concat to merge VCF
#   of non synonymous SNPs
#
#   Input:
#       - {samples}_non_synonymous.vcf
#
#   Output:
#       - all_non_synonymous.vcf
#
#   Parameters:
#       No fancy parameters here...

rule merge_vcf:
    input:
        vcfs = expand(config["VAR"] + "{samples}/{samples}_non_synonymous.vcf", samples = config["samples"])
    output:
        non_synos = config["VAR"] + "Non_synonymous_SNP.vcf"
    params:
        bind = config["BIND"],
        cont = config["CONT"]
    shell:
        """
        echo '##fileformat=VCFv4.1\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT'\
        > {output.non_synos}
        singularity exec -B {params.bind} {params.cont} \
        awk '!/#/' \
        {input.vcfs} >> {output.non_synos}
        """
