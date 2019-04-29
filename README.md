[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Documentation Status](https://readthedocs.org/projects/allmine/badge/?version=latest)](https://allmine.readthedocs.io/en/latest/?badge=latest)      [![Github All Releases](https://github.com/tbersez/Allmine)]()

# AllMine, a flexible pipeline for Allele Mining
### What is AllMine ?
AllMine is flexible and parallel pipeline for allele mining. You can use AllMine to discover *de novo* Single Nucleotide Polymorphism (**SNPs**) onto next generation sequencing data of various types (RNAseq, WGS, RRGS *ect.*). AllMine working principle is to compare your sequencing data to a annotated reference sequence to call SNPs. By defining target regions onto your reference sequence, your can extract specific SNPs and discover allele polymorphism. 
### Depandancies
* Snakemake >= V4.8.0
* Python3 >= V3.5.3
    * including modules csv and yaml
* Slurm >= V0.4.3
* Singularity >= V2.5.1
### Documentation and start guide
https://allmine.readthedocs.io/en/latest/?badge=latest
