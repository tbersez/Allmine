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

                            # Allmine, a flexible pipeline for allele mining.
                            # Develloped by Thomas Bersez (2019)
                            # INRA-GAFL / Paris Saclay university
                            # contact: thomasbersez@gmail.com


configfile: "config.yaml"
cwd = os.getcwd() + "/"
# modules loading...
include : cwd + "modules/" + config["FASTP"]
include : cwd + "modules/bwa_index_building.py"
include : cwd + "modules/" + config["BWA"]
include : cwd + "modules/varscan.py"
include : cwd + "modules/varscan_filtering.py"
#############################
rule all:
    input:
        expand(config["VAR"] + "{samples}_varscan_filtered.tab", samples = config["samples"])
