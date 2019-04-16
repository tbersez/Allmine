# Annovar format converter
#
#   This module convert vcf4 format files from Varscan
#   to annovar input format
#
#   Input:
#       - sample(vcf4)_varscan.vcf
#
#   Output:
#       - sample_varscan.avinput
#
#   Parameters:
#       None...

rule vcf4_to_avinput:
    input:
        vcf = config["VAR"] + "{samples}/{samples}_varscan.vcf"
    output:
        annovar = config["VAR"] + "{samples}/{samples}_varscan.avinput"
    params:
        bind = config["BIND"],
        cont = config["CONT"]
    message: "Converting {input.vcf} to avinput format \n"
    shell:
        """
        singularity exec -B {params.bind} {params.cont} \
        /opt/annovar/convert2annovar.pl \
        -format vcf4 \
        -includeinfo \
        -outfile {output.annovar} \
        {input.vcf}
        """
