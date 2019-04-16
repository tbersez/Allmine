#!/usr/bin/env python3
#
# This script  build the Annovar database
# for your reference genome and annotation
# config.yaml use as input

import yaml
import sys
import os

print("\33[93mAnnovar data base building may take a while.\33[0m\n")
config = yaml.safe_load(open('config.yaml'))
os.system('rm -rf avdb/')
os.system('mkdir avdb')
if config["ANO"].endswith(".gff"):
    os.system('singularity exec -B '+ config["BIND"] + ' ' + config["CONT"] + \
    " gff3ToGenePred " + config["REF"] + config["ANO"] + " avdb/AV_refGene.txt")
elif config["ANO"].endswith(".gtf"):
    os.system('singularity exec -B '+ config["BIND"] + ' ' + config["CONT"] + \
    " gtfToGenePred " + config["REF"] + config["ANO"] + " avdb/AV_refGene.txt")
os.system('singularity exec -B '+ config["BIND"] + ' ' + config["CONT"] + \
" perl /opt/annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile " + \
config["REF"] + config["GENOME"] + " avdb/AV_refGene.txt --out avdb/AV_refGeneMrna.fa")
print('\n\n\33[32mDone! Database builded successfully.\33[0m\n')
