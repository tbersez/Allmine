#!/usr/bin/python3

#   This script generates mutated sequences
#   from non synonymous SNPs found by AllMine

import sys
import yaml
import csv
import os

# reading config file (output of cvs_to_yaml.py)
config = yaml.safe_load(open('config.yaml'))
# reading non synonymous SNPs file
non_syno_snps = config["VAR"] + "Non_synonymous_SNP.vcf"
# reading bed file
with open (config["REGIONS"], 'r') as f:
    bed = [row for row in csv.reader(f,delimiter=' ')]
    # format convertion for Samtools
    regions = [bed[reg][0] + ":" + bed[reg][1]\
    + "-" + bed[reg][2] for reg in range(len(bed))]

getseq = "singularity exec -B " + config["BIND"] + " " + config["CONT"] + \
" samtools faidx " + config["REF"] + config["GENOME"] + " " + regions[0]
gene_DNA = os.popen(getseq).read()
gene_DNA = gene_DNA.split('\n')[1:]
gene_DNA = "".join(gene_DNA)
gene_DNA = gene_DNA.upper()
print(gene_DNA)
start_idx = gene_DNA.find("ATG")
for i in range(start_idx + 3, len(gene_DNA)):
    if gene_DNA[i:i+3] in ["TAG", "TAA", "TGA"] :
        print(gene_DNA[i:i+3])
        print(i)
        i += 2




sys.exit()
