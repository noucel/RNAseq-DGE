# RNAseq DGE Multiple Sample Comparision

## Description

These scripts provide a straight forward workflow to analyse RNA sequenced data from multiple datasets and compares the differential genetic expression between sample groups. RNAseq datasets can be pulled from online repositories such as NCBI Gene Expression Omnibus, Sequence Read Archive, and others.

Fasta files are to be downloaded and processed. Once this raw data is obtained, the scripts guide through several steps to perform quality control (QC) checks, for ensuring the consistency and reliability of the data. As the first step, each sequence read is assigned to its targeted gene on the genome, as a reference to knowing the specific gene that has been amplified in the experimental work. Secondly, all reads are counted genewise. Then the differential gene expression analysis is performed to observe differences between sample groups.

As a more detailed report, the resulting genes can be classified by their Gene Ontology (GO) enrichment. This is a method for investigating sets of genes using the GO system of classification, in which genes are assigned to a particular set of terms for three major domains considering cellular component, biological process, and molecular function.


## Prerequisites

Install the following R packages from Bioconductor limma, edgeR, Glimma, org.Mm.eg.db, GOstats, gplots, RColorBrewer, NMF, BiasedUrn, fastqrc, Rsubread.
Also, intall from the CRAN package repository the R packages tidyverse and DT. 


## Installation

Open the RNAseq-DGE workflow scripts.


## Workflow Overview

These work conducts a differential gene expression (DGE) analysis using RNA-seq. It starts by obtaining the raw RNA-seq data, preferable in FASTQ format, from a sequencing platform or public repositories like NCBI SRA or ENA.

Process a QC to assess and filter low-quality reads. Trim data to remove adapters and low-quality bases.

Align and map reads to the reference genome or transcriptome for which the RNA guide library has been designed.

Quantify and count the number of reads mapping to each gene or transcript. Normalise and ajust for variations in sequencing depth and library size.

Conduct the DGE analysis to identify genes or transcripts that are differentially expressed between experimental conditions (e.g., treatment vs. control). The package alograithm variations used are those included in the edgeR package.

Perform statistical tests to assess differential expression.

Correct for multiple testing by adjusting p-values for multiple comparisons using methods like Benjamini-Hochberg correction.

Annotate differentially expressed genes with functional information.

Visualise data by creating plots (e.g., volcano plots, heatmaps) to interprete the results and observe their distribution or association.

Interpret the biological significance of the differentially expressed genes and their potential roles in the studied process or condition.

Summarize results in a report or publication, including methods, results, and figures.


## Data

Data sets are best processed from raw sources of renowned data archive repositories. Contat the persons responsible for conducting the experimental work if necessary clarifing the method in which the samples were processed. Formats should be in fasta and fastq files. To obtain such dataformats pull datasets from open well-known RNA-seq repositories like those noted in the references below.


## Results

The results of this RNAseq DGE anlaysis need to be assessed in the context of each specific experimental purpose. Discussing the results with professional biologists or medical scientists is highly recommended for best interpretations.


## Contributing

Contribute to this project, report possible issues for the improvment of this workflow. Suggesting improvements is encouraged.


## License
GNU General Public License (GPL). The GPL is designed to ensure that software remains open source and freely available for anyone to use, modify, and distribute while also requiring any derivative works to be licensed under the same terms.


## Contact
Git Hub author Noucel.


## Acknowledgments

Thanks to the advice recived by colleagues and The Unviersity of Melbourne Bioinformatics gruop.


## References

NCBI Gene Expression Omnibus (GEO)
Website: https://www.ncbi.nlm.nih.gov/geo/
GEO is a comprehensive repository hosted by the National Center for Biotechnology Information (NCBI). It contains a vast collection of gene expression data, including RNA-seq datasets.

Sequence Read Archive (SRA)
Website: https://www.ncbi.nlm.nih.gov/sra
SRA is also part of the NCBI and focuses on storing raw and processed data sequences, including RNA-seq data.
