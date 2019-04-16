        ##########################################################################################
        #      ___           ___       ___       ___                       ___           ___     #
        #     /\  \         /\__\     /\__\     /\__\          ___        /\__\         /\  \    #
        #    /::\  \       /:/  /    /:/  /    /::|  |        /\  \      /::|  |       /::\  \   #
        #   /:/\:\  \     /:/  /    /:/  /    /:|:|  |        \:\  \    /:|:|  |      /:/\:\  \  #
        #  /::\~\:\  \   /:/  /    /:/  /    /:/|:|__|__      /::\__\  /:/|:|  |__   /::\~\:\  \ #
        # /:/\:\ \:\__\ /:/__/    /:/__/    /:/ |::::\__\  __/:/\/__/ /:/ |:| /\__\ /:/\:\ \:\__\#
        # \/__\:\/:/  / \:\  \    \:\  \    \/__/~~/:/  / /\/:/  /    \/__|:|/:/  / \:\~\:\ \/__/#
        #      \::/  /   \:\  \    \:\  \         /:/  /  \::/__/         |:/:/  /   \:\ \:\__\  #
        #      /:/  /     \:\  \    \:\  \       /:/  /    \:\__\         |::/  /     \:\ \/__/  #
        #     /:/  /       \:\__\    \:\__\     /:/  /      \/__/         /:/  /       \:\__\    #
        #     \/__/         \/__/     \/__/     \/__/                     \/__/         \/__/    #
        ##########################################################################################


                            # AllMine, a flexible pipeline for Allele Mining #
                            # This software is under the MIT License         #
                            # Copyright Thomas Bersez  2019                  #
                            # INRA-GAFL / Paris Saclay university            #
                            # contact: thomasbersez@gmail.com                #


#   This Snakefile script controls the loading and execution of the different
#   modules/<name>.py. It is use to start the AllMine pipeline.
#   The "config.yaml" file, generated by csv_to_yaml.py, must be present in the
#   working directory.
#   Several directories containing, processed reads, bam and putative variants
#   will be created in the working directory by the pipeline.
#   Make sure that you have the needed rights, space and resources
#   where you run AllMine!


configfile: "config.yaml"
cwd = os.getcwd() + "/"

# modules loading...
include : cwd + "modules/" + config["QC"]
include : cwd + "modules/" + config["FASTP"]
include : cwd + "modules/" + config["INDEXER"]
include : cwd + "modules/" + config["ALLIGNER"]
include : cwd + "modules/bam_with_bed_parse.py"
include : cwd + "modules/vcf_to_avinput.py"
include : cwd + "modules/varscan.py"
include : cwd + "modules/annovar.py"
include : cwd + "modules/make_report.py"
include : cwd + "modules/bam_index.py"
include : cwd + "modules/whatshap.py"

# target files...
rule all:
    input:
        expand(config["VAR"] + "{samples}/{samples}_varscan_phased.vcf", samples = config["samples"]),
        'Non_synonymous_variants_summary.tab'
