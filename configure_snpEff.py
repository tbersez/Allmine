#!/usr/bin/env python3
#
#   This file completes the snpEff.config file
#   with proper arguments in order to make it
#   works on your data. It must be run in the
#   Snakefile directory.

import yaml
import sys
import os

# loading the config.yaml file created by csv_to_yaml.py
config = yaml.safe_load(open('config.yaml'))
os.system('rm -rf snpEff_db/')
os.system('mkdir snpEff_db')

# copy files in snpEff_db
cp_genome = 'cp ' + config['REF'] + config['GENOME'] + ' snpEff_db/sequences.fa'
os.system(cp_genome)
cp_anno = 'cp ' + config['REF'] + config['ANO'] + ' snpEff_db/genes.gff'
os.system(cp_anno)

# edit snpEff.config
os.system('echo ' + "snpEff_db.genome : snpEff_db" + ' >> snpEff.config')

# build the db
os.system('singularity exec -B /mnt/nas_eic/gafl01/home/gafl/tbersez ~/Allmine/AllMine java -jar /snpEff/snpEff.jar build -v -gff3 snpEff_db')

sys.exit()
