# Tn5 library analysis
## Introduction

In the transposon mutant libraries pooled mutant strains are marked with randomly generated unique barcode sequences that allow their identification. These sequences are introduced into the genome using Tn5 or mariner transposones. Due to these sequences, we are able to estimate the abundance of each mutant strain grown together before and after cultivation via Illumina sequencing. The ratio of strain abundancies after and before cultivation is commonly referred as strain fitness.

Transposon mutant libraries are useful tools for fitness profiling in bacteria and allow to find out which mutations suppress or stimulate bacterial growth. Here, we introduce such a library for *Escherichia coli* K-12 (substrain W3110 RpoS<sup>+</sup>). The library was created and annotated as described in the following publication:

*Wetmore KM, Price MN, Waters RJ, Lamson JS, He J, Hoover CA, et al. Rapid quantification of mutant fitness in diverse bacteria by sequencing randomly bar-coded transposons. mBio. 2015;6(3):e00306-15.
https://journals.asm.org/doi/10.1128/mbio.00306-15*

Please refer to this publication for all the methodological details.

We wrote a set of scripts to use when mapping the library and using it in the actual experiments. To start annotating, one needs sequencing results and a good genome annotation of an organism they work with.

## Set of scripts to analyse Tn5 transposon libraries

1. In the part one of the library annotation, the script (*tnseq_phase1.py*) looks for the barcode flanking sequences in each read of the FASTQ file, and this is able to extract the barcodes and perform their quality check. After passing the QC, the reads containing barcodes are stored in the separate FASTA file, and the information about the barcodes is also collected in the csv file.
2. The selected reads are then mapped to the reference genome using *ncbi-blast* toolkit (*BLASTN_command.sh*).
3. Finally, the positions of each barcode are determined, and a comprehensive annotation is created (*tnseq_phase2.py*).

| seq | direction | position | gene |
| --- | --- | --- | --- |
| ATTAGAGTGGTAGTACTGGC | 1 | 248 | thrL |
| TGGTTGTGGGCAGGCGGGGT | 1 | 16872 | mokC |

When the library is annotated, one is ready to analyze the actual experiments using *barcode.counting.py* to count the abundance of each barcode in the sequencing data, and then *TnSeqW3110.R* to calculate the actual fitness data.

