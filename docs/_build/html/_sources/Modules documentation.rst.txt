Modules documentation
=====================

AllMine is developed in a **modular fashion**. Current AllMine modules are
documented here. Some of them include **tweakable parameters** ! Only key steps
modules are documented here.

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
  * ``--max_len1 350`` : Maximum lenght of submited reads. Longer reads will be
    trimmed from 3' end to the maximum size.

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
To perfom a two pass mapping strategy

BWA mem; DNAseq read mapping
----------------------------

If you have submited DNAseq data, AllMine will map your reads using
the ``bwa mem`` command from `BWA <http://bio-bwa.sourceforge.net/>`_ v0.7.17.
Two modules, **bwa_pe** and **bwa_se** have been implemented to handle
paired end and single end reads respectively. This mapping algorithm handle
reads from 70bp to 1Mbp.

Parameters :

  * ``-w`` : Set to 100 bp. Band width. Essentially, gaps longer than INT will
    not be found.
  * ``-d`` : Set to 100. Z-dropoff. Avoids unnecessary extension, but also
    reduces poor alignments inside a long good alignment.
  * ``-r`` : Set to 1.5. **Larger value yields fewer seeds, which leads to faster
    alignment speed but lower accuracy**.
  * ``-B`` : Set to 4. Mismatch penalty.
  * ``-O`` : Set to 6. Gap opening penalty.
  * ``-E`` : Set to 1. Gap extention penalty.
  * ``-L`` : Set to 5. When performing SW extension, BWA-MEM keeps track of the
    best score reaching the end of query.

Note : In most of the cases BWA parameters are well suited for balance between
accuracy and computational cost.

STAR; RNAseq reads mapping
--------------------------

If you have submited RNAseq data, AllMine will map your reads using ``STAR``
command from `STAR <https://github.com/alexdobin/STAR>`_ v2.7.0f. AllMine apply
a two pass mapping strategy. Two both handle paired and single end data four
modules where developped **star_pe_FP**, **star_pe_SP**, **star_pe_FP** and
**star_pe_SP**.

Parameters :

  * ``--scoreGap`` : Set to 0. Gap penalty.
  * ``--scoreGapNoncan`` : Set to -8. Non-canonical junction penalty.
  * ``--scoreGapGCAG`` : Set to -4. GC/AG and CT/GC junction penalty.
  * ``--scoreGapATAC`` : Set to -8. AT/AC and GT/AT junction penalty.
  * ``--scoreGenomicLengthLog2scale`` : Set to -0.25. Extra score
    logarithmically scaled with genomic length of the alignment :
    scoreGenomicLengthLog2scale*log2(genomicLength).
  * ``--scoreDelOpen`` : Set to -2. Deletion open penalty.
  * ``--scoreDelBase`` : Set to -2. Deletion extension penalty per base
    (in addition to scoreDelOpen).
  * ``--scoreInsOpen`` : Set to -2. Insertion open penalty.
  * ``--scoreInsBase`` : Set to -2. Insertion extension penalty per base
    (in addition to scoreInsOpen).
  * ``--scoreStitchSJshift`` : Set to 1. Maximum score reduction while searching
    for SJ boundaries in the stitching step.
  * ``--runThreadN`` : Set to 10. Number of threads used per jobs.

Note : STAR include numerous parameters, please read the manual for more
informations.

Varscan; SNP calling
--------------------

AllMine use `Varscan <http://varscan.sourceforge.net/>`_ v2.4.3 to perform SNP
calling on aligned reads.

Parameters :

  * ``--p-value`` : Set to 0.99. We did not choosed to use p value as filter for
    SNP calling. However you can tune this parameter if wanted.
  * ``--min-coverage`` : Set to 8. Minimal deep at SNP loci for calling.
  * ``--min-var-freq`` : Set to 0.15. Minimal variant frequency at SNP loci for
    calling.
  * ``--min-avg-qual`` : Set to 20. Minimal sequencing quality at SNP loci for
    calling.

Note : The tunning of Varscan parameters is a trade-off between sensitivity and
specificity. It should be keeped in mind when analysing your results!

WhatsHap; SNP Phasing
---------------------

AllMine use `WhatsHap <https://bitbucket.org/whatshap/whatshap>`_ v0.18 to
perform phasing of called SNPs.

Note : 330 Bp long paired end reads at 15X depth are required for haplotypes
chunks of 300 Kb in average.

ANNOVAR; SNP annotation
-----------------------

AllMine use `SnpEff <http://annovar.openbioinformatics.org/en/latest/>`_
v2018-04-16 00:43:31 to annotate called SNPs. A gene base annotation is done by
default.
