Quick Start guide
=================

First check that all requiered :ref:`to-dep` are installed. If you wish to use
AllMine on a cluster, contact your admnistrator for instalations.

Getting AllMine code
--------------------

AllMine code is GitHub hosted. To get the source code using git CLI, use :

.. code-block:: bash

    git clone https://github.com/tbersez/Allmine.git

Then, from within AllMine directory, build **AllMine container** using
Singularity (admin rights needed) :

.. code-block:: bash

    sudo singularity build AllMine singularity/AllMine_recipe.sr

Here, AllMine is ready to run.

Input files
-----------

To perform allele mining, AllMine needs :

 * Sequencing data (DNA or RNA)
 * A reference genome and annotation
 * A bed file with regions (or genes) of interest

Sequencing data
---------------

AllMine supports **RNA and DNA** sequencing, **paired or single end**, with a
maximum read lenght of **350bp** (you can submit longer reads but they will be
trimmed to the maximum size). Sequencing data must be in **fastq format** and
may be gzip compressed.

Sample sheet
------------

Configuration maker for AllMine use a **csv** describing samples file as input.
The csv file must have the flowing headers :

.. code-block:: bash

        filename R1_ext R2_ext platform date(mm/dd/yy)

Here is an example:

.. code-block:: bash

        "filename","R1_ext","R2_ext","platform","date(mm/dd/yy)"
        "SRR1538456","_1.fastq.gz","_2.fastq.gz","illumina","03/03/03"
        "SRR1538457","_1.fastq.gz","_2.fastq.gz","illumina","04/03/03"
        "SRR1538484","_1.fastq.gz","_2.fastq.gz","illumina","05/03/03"

You can use a spread sheet editor to create this csv file !

Reference genome and annotation
-------------------------------

Reference genome must be provided in one file, in **fasta format**.
Annotation can be provided in **gff** or **gtf** format (recommended).
When possible, we advice you to download the reference sequence and annotation
from curated sources, such as `Ensembl <http://ensemblgenomes.org/>`_.

Bed file
--------

Regions of interest must be specified using a **bed** file, here is an example :

.. code-block:: bash

        NC_035163.1 25395963 25398308
        NC_035168.1 31042453 31045884
        NC_035169.1 25941228 25944616
        NC_035175.1 3317633 3320503
        NC_035177.1 20184932 20187543

Be sure that the first column (contigs) match with the ones used in your
reference genome.

Making your work space ready
----------------------------

Place your reference genome and annotation in a common folder. That one must
only contain those both files. Place your bed file with regions of interest in
an other folder. **AllMine outputs are created where your start the analysis.**
Make sure that you have enought space to store all outputs !

Configuration of an AllMine run
-------------------------------

To configure an AllMine run use :

.. code-block:: bash

          ./csv_to_yaml.py path/to/sample_sheet.csv

Answer the questions the script is asking you to configure your run.
Note : the bind path is the path from the **root to your home folder.**

Once done run :

.. code-block:: bash

          ./annovar_makebd.py

This script build the annotation database of Annovar. It need to done once for
each new genome used.

Running AllMine
---------------

We recommend first to do a dry run using the following command.
**CORE_NUMBER** must be replaced by the number of cores you wish to use.

.. code-block:: bash

          snakemake -j CORE_NUMBER \
          --cluster-config slurm_config.json \
          --cluster "sbatch" -n

Check the output to ensure that your run is properly configured.
If not, return to configuration step to correct errors. If yes, run :

.. code-block:: bash

          snakemake -j CORE_NUMBER \
          --cluster-config slurm_config.json \
          --cluster "sbatch"

AllMine is now running. Depending on how much data you have submited and
CORE_NUMBER, the analysis may take from few hours to a few days.
