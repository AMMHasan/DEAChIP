# DEAChIP

**D**ifferential **E**nrichment **A**nalysis of **ChIP**-seq data

This script describes the workflow for comparative ChIP-seq data (from Illumina sequence platform) analysis, to accompany the manuscript "Genomic Analysis of DNA Double-Strand Break Repair in *Escherichia coli*" (January 2018. Methods in Enzymology; DOI:10.1016/bs.mie.2018.09.001)

The analysis can be done manually with a combination of library size normalisation (suggested by *Simon Anders & Huber, 2010*) and further normalisation of IP data by input data (where available).

In addition, R/Bioconductor package DESeq is used here for the comparative analysis of IP data (with biological replicates) between different conditions or different strains of *Escherichia. coli*. In that case, at least three biological replicates of each strain or condition are needed for a statistically significant fold-change measurement of DNA enrichment at a defined region in the genome, otherwise, the same method can be used for exploratory data analysis without any statistical significance. In the latter case, a minimum of two independent biological repeats are recommended to confirm qualitative reproducibility.

