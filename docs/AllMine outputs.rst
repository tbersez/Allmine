AllMine outputs
===============

AllMine outputs are detailled here. Some of them (such a trimmed reads)
are "non-final" outputs and so, should not be conserved after the completion
of the run. However if, for any reasons you may need them, you can gather them
from the outputs.

Quality control reports and trimmed reads
-----------------------------------------
During prepossessing, adapters, linker, low quality bases and low complexity
regions are trimmed for the raw reads. Trimmed reads are stored in the
**trimmed** directory. Thoses reads are then passed to quality checking.
Quality control steps are performed by AllMine in order to ensure the quality of
the sequencing data submitted. Reports names use the format <sample_name>.html
and can be found in the **QC_post_preproc** directory.
You can inspect them using your favorite web browser. As the sequencing data
quality influence allele mining reliability, quality control report always
should be carefully inspected and taken in account during the results
interpretation.

Mapped reads and de novo splicing junctions
-------------------------------------------

Two different strategies of mapping are supported by AllMine. For DNA sequencing
inputs, AllMine uses BWA (Burrows-Wheeler Aligner). For RNA sequencing input,
STAR (Spliced Transcripts Alignment to a Reference) is used with a two pass
strategy. In both cases, mapped and sorted reads are stored in the **mapped**
directory. Both full alignments files and parsed around specified regions are
conserved, allowing new parsing around other regions without executing the
mapping once again. In the case of RNA sequencing data, the sub-directory
**STAR_SJ** is created. It contains the de novo splicing junctions discovered by
STAR during the first pass of the mapping.

Putative SNPs
-------------

Putative variants called by Varscan can be found in the **variant** directory.
Each subfolder correspond to one sample. Annotated, phased and raw variants are
displayed.

Non synonymous variant summary
------------------------------

This is the main AllMine output. **Non_synonymous_variant_summary.tab** is a
tabular file displaying all non synonymous SNPs found by AllMine. Base and amino
acid reference and variation are indicated as well as genotype (het or hom),
location in the gene (exon1, exon2 ...). If a variant is found in more than one
sample, all concerned samples are indicated in the SAMPLE(s) column.

Run Report
----------

Allmine build an R markdown report to sum up most of the information about
the run. This report also include coverage plots for regions of interest.
