Modules documentation
=====================

AllMine is developed in a **modular fashion**. Current AllMine modules are
documented here. Some of them include tweakable parameters !

fastp; reads presprocessing
---------------------------

AllMine uses `fastp <https://github.com/OpenGene/fastp>`_ v0.20.0 to
perform presprocessing on submited reads. Two modules, **fastp_pe** and
**fastp_se** have been implemented to handle paired end and single end reads
respectively.

Parameters :

  * ``--correction`` : Enable paired end overlaping regions correction.
  * ``--cut_mean_quality`` : Set to 20. Threshold for trimming in **both 5' and 3'**.
    You can increase this number for a more stringeant quality trimming.
  * ``--cut_window_size`` : Size of the sliding window for sliding window trimming.
    Increasing this number will relax the trimming. Set to 1 for per base trimming
    (not recommended).
  * ``--complexity_threshold`` : Set to 30. Complexity is defined as **P(base[i] !=
    base[i+1])**. Reads bellow the threshold are discarded. Usefull to filter polyA
    tails in RNA seq data !
  * ``-w`` : Threads used.

Please refer to the fastp manual for more informations.

FastQC; quality control
-----------------------

AllMine uses `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
v0.11.8 to perform quality control on submited reads. Two modules, **fastqc_pe**
and **fastqc_se** have been implemented to handle paired end and single end reads
respectively. Numerous indices and statistics are computed by FastQC. You may
inspect them to ensure your data quality.

BWA index; genome indexing
--------------------------

If you have submited DNAseq data, AllMine will index your reference genome using
the ``bwa index`` command from `BWA <http://bio-bwa.sourceforge.net/>`_ v0.7.17.

Parameters :

  * ``-a is`` : Algorithm used to index the genome. **IS** does not work for
    large genome (size > 3Gbp). Switch to ``-a bwtsw`` if needed.


STAR index; genome indexing
---------------------------

If you have submited RNAseq data, AllMine will index your reference genome using
the ``STAR --runMode genomeGenerate`` command from
`STAR <https://github.com/alexdobin/STAR>`_ v2.7.0f.
