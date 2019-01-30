#!/usr/bin/env python3

#   This script allow easy generation of .yaml configuration files
#   via your bash promt. Please note that all directories paths must
#   specified folowing the <directory_name>/ format, and so need to
#   be in the working directory.


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
            yaml.write("TECH: RNAseq \n")
        elif(tech == 'r'):
            yaml.write("TECH: DNAseq \n")
        if(mode == 'p'):
            yaml.write("MODE: paired \n")
            yaml.write("FASTP: fastp_pe.py \n")
            yaml.write("BWA : bwa_pe.py \n")
        elif(mode == s):
            yaml.write("MODE: single \n")
            yaml.write("FASTP: fastp_se.py \n")
            yaml.write("BWA: bwa_se.py \n")
        yaml.write("RAW: " + cwd + reads_dir + "\n")
        yaml.write("TRIMMED: " + cwd + "trimmed/ \n")
        yaml.write("REF: " + cwd + genome_dir + "\n")
        yaml.write("GENOME: " + os.listdir(genome_dir)[0] + "\n") # only one reference may be provided!
        yaml.write("MAP: " + cwd + "mapped/ \n")
        yaml.write("VAR: " + cwd + "variant/ \n")
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
