#!/usr/bin/python3

#   This script generates a summary of all non synonymous
#   SNPs found by AllMine in the regions of interest.

import sys
import yaml

# reading config file (output of cvs_to_yaml.py)
config = yaml.safe_load(open('config.yaml'))
# list of variants
de_novo_var = []
# associated samples (max size 1000)
is_in_sample = [None] * 1000
samples = list((config["samples"]).keys())
for sample in samples:
    # path to vcf with non synonymous SNPs
    path = config['VAR'] + sample + "/" + sample + "_non_synonymous.vcf"
    with open (path, 'r') as vcf:
        for lines in vcf.readlines():
            # skip comment lines
            if lines[0] == "#":
                continue
            # parsing vcf lines
            line = (lines.strip().split("\t"))
            chr = line[0]
            base = line[1]
            ref = line[3]
            alt = line [4]
            eff = line[7].strip().split("|")[1]
            gene = line[7].strip().split("|")[18]
            variant = (chr, base, ref, alt, gene)
            if variant not in de_novo_var:
                # new variant
                de_novo_var.append(variant)
                id = de_novo_var.index(variant)
                is_in_sample[id] = sample
            else:
                id = de_novo_var.index(variant)
                if is_in_sample[id] == None:
                    # create entry for new variant
                    is_in_sample[id] = sample
                else:
                    # build sample list for variant of index 'id'
                    is_in_sample[id] = is_in_sample[id] + (", " + sample)
with open('Non_synonymous_variants_summary.tab', 'w') as out:
    id = 0
    out.write("chr\tposition\tref\talt\tgene\tsamples\n")
    for snp in de_novo_var:
        for i in snp :
            out.write(i + '\t')
        out.write(is_in_sample[id])
        out.write('\n')
        id += 1
sys.exit()
