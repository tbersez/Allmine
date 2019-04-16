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
    path = config['VAR'] + sample + "/" + sample + "_varscan.avinput.exonic_variant_function"
    with open (path, 'r') as vcf:
        for lines in vcf.readlines():
            # skip comment lines
            if lines[0] == "#":
                continue
            # parsing vcf lines
            line = (lines.strip().split("\t"))
            if(line[1] == 'nonsynonymous SNV'):
                chr = line[3]
                base = line[4]
                ref = line[6]
                alt = line[7]
                aa_change = line[2].strip().split(":")[4].split(",")[0].split(".")[1]
                genid = line[2].strip().split(":")[0]
                genotype = line[8]
                exon = line[2].strip().split(":")[2]
                variant = (chr, base, ref, alt, aa_change, exon, genid, genotype)
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
                        is_in_sample[id] = is_in_sample[id] + ("," + sample)
with open('Non_synonymous_variants_summary.tab', 'w') as out:
    id = 0
    out.write("CHR\tBASE\tREF\tALT\tAA_CHANGE\tEXON\tGENE\tGENOTYPE\tSAMPLE(s)\n")
    for snp in de_novo_var:
        for i in snp :
            out.write(i + '\t')
        out.write(is_in_sample[id])
        out.write('\n')
        id += 1
sys.exit()
