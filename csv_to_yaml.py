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

# TODO: add config.yaml removal before error exit

import sys
import csv
import os

tab = '    '
cwd = os.getcwd() + "/"

reads_dir = input('Path to raw read directory : ')
genome_dir = input('Path to the reference genome_dir : ')
tech = input('DNAseq or RNAseq (enter d or r) : ')
mode = input('Paired end or single end (enter p or s) : ')
threads = input('How many threads should be used : ')

with open(sys.argv[1], 'r') as file:
    reader = csv.DictReader(file, delimiter=',', quotechar='"')
    with open("config.yaml", 'w') as yaml:
        if(tech == 'd'):
            yaml.write("TECH: DNAseq \n")
            files = os.listdir(cwd + genome_dir)
            for file in files: # Only one fasta is accepted !
                if file.endswith(".fna" or ".fa" or ".fasta"):
                    yaml.write("GENOME: " + file + "\n")
                else:
                    print("No valid reference provided in " + genome_dir + " !\n")
                    sys.exit()
        elif(tech == 'r'):
            files = os.listdir(cwd + genome_dir)
            for file in files: # Only one fasta and one gff are accepted !
                if file.endswith(".fna" or ".fa" or ".fasta"):
                    yaml.write("GENOME: " + file + "\n")
                elif file.endswith(".gff"):
                    yaml.write("ANO: " + file + "\n")
                else:
                    print("Non valid files in " + genome_dir + " !\n")
                    sys.exit()
            yaml.write("TECH: RNAseq \n")
        if(mode == 'p'):
            yaml.write("MODE: paired \n")
            yaml.write("FASTP: fastp_pe.py \n")
            yaml.write("BWA : bwa_pe.py \n")
        elif(mode == 's'):
            yaml.write("MODE: single \n")
            yaml.write("FASTP: fastp_se.py \n")
            yaml.write("BWA: bwa_se.py \n")
        yaml.write("REF: " + cwd + genome_dir + "\n")
        yaml.write("RAW: " + cwd + reads_dir + "\n")
        yaml.write("TRIMMED: " + cwd + "trimmed/ \n")
        yaml.write("MAP: " + cwd + "mapped/ \n")
        yaml.write("VAR: " + cwd + "variant/ \n")
        yaml.write("LOG: " + cwd + "log/ \n")
        yaml.write("THREADS: " + threads + "\n\n")
        yaml.write("samples:\n")
        for row in reader:
            yaml.write(tab + row["filename"] + ": \n")
            yaml.write(2*tab + "name: " + row["filename"] + "\n")
            yaml.write(2*tab + "R1: " + cwd + reads_dir + row["filename"] + row["R1_ext"] + "\n")
            if(mode == 'p'):
                yaml.write(2*tab + "R2: " + cwd + reads_dir + row["filename"] + row["R2_ext"] + "\n")
            yaml.write(2*tab + "platform: " + row["platform"] + "\n")
            yaml.write(2*tab + "date: " + row["date(mm/dd/yy)"] + "\n")
sys.exit()

# TODO: add config.yaml remove
