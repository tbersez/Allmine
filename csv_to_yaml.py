#!/usr/bin/env python3

#   This script allow easy generation of .yaml configuration files
#   via your bash promt. Please note that all directories paths must
#   specified folowing the <directory_name>/ format, and so need to
#   be in the working directory. Please note that one and only one
#   reference sequence (and gff for RNAseq) are accepted at a time !
#   Reference must have the extention ".fna", ".fa" or ".fasta".
#   Annotation must have the extention ".gff".
#   If one, or more, of this conditions are not filled, this script
#   will exit with an error message.

import sys
import csv
import os
from os.path import expanduser

tab = '    '
cwd = os.getcwd() + "/"
bind = expanduser("~")

# AllMine !
print ("""\33[1m\33[36m
        ##################################################
        # AllMine, a flexible pipeline for Allele Mining #
        # This software is under the MIT License         #
        # Copyright Thomas Bersez  2019                  #
        # INRA-GAFL / Paris Saclay university            #
        # contact: thomasbersez@gmail.com                #
        ##################################################
        \33[0m
       """)


# User interface
print("\33[93mPlease answer to the following questions to configure your AllMine run.")
print("Paths must end with a /. Press enter to confirm.\33[0m\n\n")
reads_dir = input('Path to raw read directory : ')
genome_dir = input('Path to the reference genome_dir : ')
regions = input('Path to bed file with regions of interest : ')
tech = input('DNAseq or RNAseq data (enter d or r) : ')
mode = input('Paired end or single end reads (enter p or s) : ')
bind = input('Bind path point to your home directory (absolute path) : ')
container_path = input('Absolute path to AllMine container : ')

# Parse the samples sheet in csv format
with open(sys.argv[1], 'r') as file:
    reader = csv.DictReader(file, delimiter=',', quotechar='"')
    with open("config.yaml", 'w') as yaml:

        # DNA mode
        if(tech == 'd'):
            yaml.write("TECH: DNAseq \n")
            files = os.listdir(cwd + genome_dir)
            for file in files: # Only one fasta and one gff are accepted !
                if file.endswith(".fna" or ".fa" or ".fasta"):
                    yaml.write("GENOME: " + file + "\n")
                elif file.endswith(".gff" or ".gtf"):
                    yaml.write("ANO: " + file + "\n")
                else:
                    print("\33[31mNo valid reference provided in " + genome_dir + " !\n\33[0m")
                    sys.exit()

        # RNA mode
        elif(tech == 'r'):
            yaml.write("TECH: RNAseq \n")
            files = os.listdir(cwd + genome_dir)
            for file in files: # Only one fasta and one gff are accepted !
                if file.endswith(".fna" or ".fa" or ".fasta"):
                    yaml.write("GENOME: " + file + "\n")
                elif file.endswith(".gff" or ".gtf"):
                    yaml.write("ANO: " + file + "\n")
                else:
                    print("\33[31mNo valid reference provided in " + genome_dir + " !\n\33[0m")
                    sys.exit()


        # PAIRED END MODE
        if(mode == 'p'):
            yaml.write("MODE: paired \n")
            yaml.write("QC: fastqc_pe.py \n")
            yaml.write("FASTP: fastp_pe.py \n")
            if(tech == 'd'):
                yaml.write("INDEXER: bwa_index_building.py \n")
                yaml.write("ALLIGNER: bwa_pe.py \n")
            elif(tech == 'r'):
                yaml.write("INDEXER: star_index_building.py \n")
                yaml.write("ALLIGNER: star_pe_FP.py \n")
        # SINGLE END MODE
        elif(mode == 's'):
            yaml.write("MODE: single \n")
            yaml.write("QC: fastqc_se.py \n")
            yaml.write("FASTP: fastp_se.py \n")
            if(tech == 'd'):
                yaml.write("INDEXER: bwa_index_building.py \n")
                yaml.write("ALLIGNER: bwa_se.py \n")
            elif(tech == 'r'):
                yaml.write("INDEXER: star_index_building.py \n")
                yaml.write("ALLIGNER: star_se_FP.py \n")

        # Paths
        yaml.write("BIND: " + bind + "\n")
        yaml.write("REF: " + cwd + genome_dir + "\n")
        yaml.write("RAW: " + cwd + reads_dir + "\n")
        yaml.write("TRIMMED: " + cwd + "trimmed/ \n")
        yaml.write("MAP: " + cwd + "mapped/ \n")
        yaml.write("VAR: " + cwd + "variant/ \n")
        yaml.write("REGIONS: " + cwd + regions + "\n")
        yaml.write("CONT: " + cwd + container_path + "\n")

        # SAMPLE objects creation
        yaml.write("samples:\n")
        for row in reader:
            yaml.write(tab + row["filename"] + ": \n")
            yaml.write(2*tab + "name: " + row["filename"] + "\n")
            yaml.write(2*tab + "R1: " + cwd + reads_dir + row["filename"] + row["R1_ext"] + "\n")
            if(mode == 'p'):
                yaml.write(2*tab + "R2: " + cwd + reads_dir + row["filename"] + row["R2_ext"] + "\n")
            yaml.write(2*tab + "platform: " + row["platform"] + "\n")
            yaml.write(2*tab + "date: " + row["date(mm/dd/yy)"] + "\n")
print('\n\n\33[32mDone. Configuration file \"config.yaml\" created in the current directory.\33[0m\n')
sys.exit()
