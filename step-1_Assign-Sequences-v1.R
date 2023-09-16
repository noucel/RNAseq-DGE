##################################################################
# The University of Melbourne
# The Peter Doherty Institute for Infection and Immunity
# Department of Microbiology and Immunology
# D. Godfrey and A. Uldrich laboratories
#
# Author: Marc Rigau
# Title: Comparison of RNA-seq between Human Cell Lines that Differential Bind to BTN-tetramers
# Date: Jan 10, 2022
# Updated: May 12, 2022
# Version: 4
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Notes

# https://monashbioinformaticsplatform.github.io/RNAseq-DE-analysis-with-R/RNAseq_DE_analysis_with_R.html
# https://www.alzheimersworkbench.ucsd.edu/EndToEndAnalysis_RNASeq.html

##################################################################
# Install the packages from Bioconductor
# if (!requireNamespace("BiocManager"))
# install.packages("BiocManager")
# BiocManager::install( c( "limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn", 'fastqcr', 'Rsubread' ) )
# or, download the installer script
#source("http://bioconductor.org/biocLite.R")

##################################################################
### Libraries
##################################################################

library( tidyverse ) # Kernel for statistics in science fields.
library( Rsubread ) # Map the reads to the reference genome (eg. the human genome “hg19”) (Y et al. 2013).

##################################################################
### Data Pre-processing
##################################################################
# Raw sequencing data are in Fastq format which is a well defined text-based format for storing nucleotide sequences with corresponding quality score.

# Define the directory of the RNAseq datasets
RNAseqFastq.dir <- "~/UniMelb/Datasets/RNAseq/Fastq-Raw"
#RNAseqFastq.dir <- "/Volumes/WolfeCreek/Marc/RNAseq/Fastq-Raw"

# List all fastq files
dir( RNAseqFastq.dir, include.dirs=T )


##################################################################
### Mapping
##################################################################

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Note
# Rsubread provides reference genome indices for the most common two organisms: human and mouse.
# If working with a different organism, build a new index using the buildindex command.
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# Mapping the sequence reads to a build a default base-space index.
# You can provide a list or single Fasta files.
# Compressed gzipped files are read.

# Name the reference genome raw fasta file
reference.genome.file.name <- "genome_GRCh38_p13.fa"

# Define the output directory for the Rsubread index
# Note: the folder have to be created beforehand for the buildindex formula 
ref.genome.dir <- "~/UniMelb/Datasets/RNAseq/GENCODE-Human_Release-40-Index"

# Define the basename for the index
human.genome <- "GRCh38.p13"


#=================================================================
# Built a hash table index for the reference genome
# This only needs to be made once
# To construct the index takes time. Expect an hour or two

# Build the index
# buildindex( basename = file.path( ref.genome.dir, human.genome ),
#             reference = file.path( ref.genome.dir, reference.genome.file.name ),
#             gappedIndex	= T
#             )
#=================================================================

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Align the reads to the reference genome in the form of a Binary Alignment Map (BAM) which is the most comprehensive raw data of genome sequencing.
# BAM files consists of lossless, compressed binary representation of the Sequence Alignment Map-files.
# BAM is the compressed binary representation of SAM, a compact and index-able representation of nucleotide sequence alignments.
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# Define the Raw data files for the forward read files
read.files.fwd <- list.files( RNAseqFastq.dir, pattern = "1.fastq.gz$" )
read.files.fwd.path <- file.path( RNAseqFastq.dir, read.files.fwd )

# Define the Raw data files for the reverse read files
read.files.rvs <- list.files( RNAseqFastq.dir, pattern = "2.fastq.gz$" )
read.files.rvs.path <- file.path( RNAseqFastq.dir, read.files.rvs )

all.equal( length( read.files.fwd ), length( read.files.rvs ) )

# Define the Mapping directory for the alignment files
#mapping.dir <- "~/UniMelb/Datasets/RNAseq/Mapping"
mapping.dir <- "/Volumes/WolfeCreek/Marc/RNAseq/Mapping"


# Define the output data files
sample.name <- substr( read.files.fwd, 1, nchar( read.files.fwd) -11 )
output.BAM.file.names <- paste0( sample.name, '.bam' )

# Define the directory path of the alignment file output
output.BAM.file.names.path <- file.path( mapping.dir, output.BAM.file.names )




align( index = file.path( ref.genome.dir, human.genome ),
       readfile1 = file.path( RNAseqFastq.dir, c('SRR13479789_GSM5025500_U266_scrbl_rep1_RNAseq_Homo_sapiens_RNA-Seq_1.fastq.gz', 'SRR13049249_GSM4905405_control_K562_RNA-seq_rep1_Homo_sapiens_RNA-Seq_1.fastq.gz') ),
       readfile2 = file.path( RNAseqFastq.dir, c('SRR13479789_GSM5025500_U266_scrbl_rep1_RNAseq_Homo_sapiens_RNA-Seq_2.fastq.gz', 'SRR13049249_GSM4905405_control_K562_RNA-seq_rep1_Homo_sapiens_RNA-Seq_2.fastq.gz') ),
       type = 'rna',
       maxMismatches = 3, # Mis-matches found in soft-clipped bases are not counted.
       output_file = file.path( mapping.dir,c( 'SRR13479789_GSM5025500_U266_scrbl_rep1_RNAseq_Homo_sapiens_RNA-Seq', 'SRR13049249_GSM4905405_control_K562_RNA-seq_rep1_Homo_sapiens_RNA-Seq.fastq.gz') ),
       nthreads = 8, # Speed up the process by running several CPUs in parallel.
       output_format = "BAM"
)



#=================================================================
# Run the align command to map the reads
align( index = file.path( ref.genome.dir, human.genome ),
       readfile1 = read.files.fwd.path[2:3],
       readfile2 = read.files.rvs.path[2:3],
       type = 'rna',
       maxMismatches = 3, # Mis-matches found in soft-clipped bases are not counted.
       output_file = output.BAM.file.names.path[2:3],
       nthreads = 8, # Speed up the process by running several CPUs in parallel.
       output_format = "BAM"
)
#=================================================================

# Return the proportion of mapped reads in the output alignment file:
# Total number of input reads
# Number of mapped reads
# Proportion of mapped reads

# Summary of the aligned/mapping read sequences
propmapped( output.BAM.file.names.path )


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Note

# Some actions can be done out of R, in a shell window:
# SAM files can be converted in BAM using samtools in the terminal window.
# BAM files save disk space on your machine, compared to SAM file format.
# To convert genome aligned-mapped files:

# for FILE in *.sam; do samtools view -Shb $FILE -o $FILE.bam; done

# Similarly, the fastq files can be compressed.
# To compress fastq files?

# gzip myfile.fq

# Visualisation of the BAM files you created can also be done via a genome browser such as IGV (http://www.broadinstitute.org/igv/) after sorting and indexing of those files.
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\






