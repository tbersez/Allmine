#!/usr/bin/env python3
# This script allow easy generation of .yaml configuration files
# used by the Snakefile

import sys
import csv

tab = '    '

reads_dir = input('Path to raw read directory : ')
genome_dir = input('Path to the reference genome_dir :')
threads = input('How many threads should be used : ')

with open(sys.argv[1], 'r') as file:
    reader = csv.DictReader(file, delimiter=',', quotechar='"')
    with open("config.yaml", 'w') as yaml:
        yaml.write("RAW: " + reads_dir + "\n")
        yaml.write("REF: " + genome_dir + "\n")
        yaml.write("THREADS: " + threads + "\n\n")
        yaml.write("samples:\n")
        for row in reader:
            yaml.write(tab + row["filename"] + ": \n")
            yaml.write(2*tab + "name: " + row["filename"] + "\n")
            yaml.write(2*tab + "tech: " + row["tech"] + "\n")
            yaml.write(2*tab + "p/s: " + row["p/s"] + "\n")
            yaml.write(2*tab + "R1: " + reads_dir + row["filename"] + row["R1_ext"] + "\n")
            if(row["p/s"] == 'p'):
                yaml.write(2*tab + "R2: " + reads_dir + row["filename"] + row["R2_ext"] + "\n")
            yaml.write(2*tab + "platform: " + row["platform"] + "\n")
            yaml.write(2*tab + "date: " + row["date(mm/dd/yy)"] + "\n")
