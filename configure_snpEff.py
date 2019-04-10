#!/usr/bin/env python3
#
#   This file completes the snpEff.config file
#   with proper arguments in order to make it
#   works on your data. It must be run in the
#   Snakefile directory.

import yaml
import sys
import os

print("\33[93mSnpEff database building may take a while.\33[0m\n\n")

# loading the config.yaml file created by csv_to_yaml.py
config = yaml.safe_load(open('config.yaml'))
os.system('rm -rf snpEff_db/')
os.system('mkdir snpEff_db')

# copy files in snpEff_db
cp_genome = 'cp ' + config['REF'] + config['GENOME'] + ' snpEff_db/sequences.fa'
os.system(cp_genome)
cp_anno = 'cp ' + config['REF'] + config['ANO'] + ' snpEff_db/genes.gff'
os.system(cp_anno)
cp_anno = 'cp ' + config['REF'] + 'GCF_001659605.1_Manihot_esculenta_v6_protein.faa snpEff_db/protein.fa'
os.system(cp_anno)

# edit snpEff.config
os.system('echo ' + "snpEff_db.genome : snpEff_db" + ' >> snpEff.config')

# build the db
os.system('singularity exec -B '+ config["BIND"] + ' ' + config["CONT"] + \
' java -jar /snpEff/snpEff.jar build -gff3 -v snpEff_db')

print('\n\n\33[32mDone! Database builded successfully.\33[0m\n')
sys.exit()
