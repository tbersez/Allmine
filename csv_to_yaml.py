#!/usr/bin/env python3

import sys
import csv

tab = '    '

reads_dir = input('Path to raw read directory (must finish with a "/" !) : ')
threads = input('How many threads should be used : ')


with open(sys.argv[1], 'r') as file:
    reader = csv.DictReader(file, delimiter=',', quotechar='"')
    with open("config.yaml", 'w') as yaml:
        yaml.write("RAW: " + reads_dir + "\n")
        yaml.write("TRIMMED: trimmed_reads/\n")
        yaml.write("THREADS: " + threads + "\n")
        yaml.write("samples:\n")
        #looping over row to create a sample object for each sample
        for row in reader:
            yaml.write(tab + row["filename"] + ": \n")
            yaml.write(2*tab + "name: " + row["filename"] + "\n")
            yaml.write(2*tab + "R1: " + reads_dir + row["filename"] + row["R1_ext"] + "\n")
            yaml.write(2*tab + "R2: " + reads_dir + row["filename"] + row["R2_ext"] + "\n")
            yaml.write(2*tab + "platform: " + row["platform"] + "\n")
            yaml.write(2*tab + "date: " + row["date(mm/dd/yy)"] + "\n")
