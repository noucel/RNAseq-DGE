##################################################################
# Institution: The University of Melbourne
# Building: The Peter Doherty Institute for Infection and Immunity
# Department: Microbiology and Immunology
# Laboratory: D. Godfrey/A. Uldrich
# Position: Postdoc Research Fellow
# Author: Marc Rigau
# Title: Quality control of FastQ files
# Date: Aor 21, 2022
# Updated: May 8, 2022
# Version: 1
##################################################################
# Notes
# Useful guideline: https://cran.r-project.org/web/packages/fastqcr/readme/README.html
##################################################################
# Install the packages from Bioconductor
# if (!requireNamespace("BiocManager"))
# install.packages("BiocManager")
# BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))
##################################################################
# Libraries
library( fastqcr ) # A high throughput sequence QC analysis tool
fastqc_install() # Required to install the FastQC Tool on Unix systems (Mac OSX and Linux)
library(dplyr)

# FastQC Package Description:
#'FASTQC' is the most widely used tool for evaluating the quality of high throughput sequencing data. It produces, for each sample, an html report and a compressed file containing the raw data.
# If you have hundreds of samples, you are not going to open up each 'HTML' page.
# You need some way of looking at these data in aggregate. 'fastqcr' Provides helper functions to easily parse, aggregate and analyze 'FastQC' reports for large numbers of samples.
# It provides a convenient solution for building a 'Multi-QC' report, as well as, a 'one-sample' report with result interpretations.


##################################################################
#### Quality Control of Fastq Data
##################################################################

# The data considered for the RNAseq part of this experiment have been pulled from GEO repository and downloaded from SRA Explorer website (https://sra-explorer.info/#).
# Raw sequencing data is in FASTQ format which is a well defined text-based format for storing both biological nucleotide sequences and their corresponding quality scores.

##################################################################
# The raw data from this study lies in the directory “~/Unimelb/Datasets/RNAseq”.

# Set the work space directory
wd.dir <- "/Volumes/Noucel/UniMelb-Sync/RNAseq"

# Define the directory path for RNA-Seq data raw Fasq files
raw.fastq.dir <- file.path( wd.dir, 'Fastq-Raw')

# List the fastq files in the raw data directory
dir( raw.fastq.dir, include.dirs = T )

# Note, the fastqc function requires a single directory pathway and does not operate with multiple directory pathways.
# A single fastqc function has to be run for each directory with fastq files

##################################################################
### Fast File Quality Control Reports
##################################################################

# For a few samples a single command can provide quality control reports.
# For each sample, FastQC performs a series of tests called analysis modules.
# These modules include:

# Basic Statistics,
# Per base sequence quality,
# Per tile sequence quality
# Per sequence quality scores,
# Per base sequence content,
# Per sequence GC content,
# Per base N content,
# Sequence Length Distribution,
# Sequence Duplication Levels,
# Overrepresented sequences,
# Adapter Content
# Kmer content

# Set the path to a new output directory
fastQC.dir <- file.path( wd.dir, 'FastQC')

# Run a FastQC tool:
fastqc(fq.dir = raw.fastq.dir, # FASTQ files directory
       qc.dir = fastQC.dir, # Results direcory
       threads = 32 # Number of threads
       )


##################################################################
### Aggregating Reports
##################################################################

# Aggragate files in the directory containing zipped FastQC reports
qc <- qc_aggregate( fastQC.dir )
head( qc )

# The table shows, for each sample, the names of tested FastQC modules, the status of the test, as well as, some general statistics including the number of reads, the length of reads, the percentage of GC content and the percentage of duplicate reads.

# Column names:
# sample: sample names
# module: fastqc modules
# status: fastqc module status for each sample
# tot.seq: total sequences (i.e.: the number of reads)
# seq.length: sequence length
# pct.gc: percentage of GC content
# pct.dup: percentage of duplicate reads

##################################################################
# Inspect modules that failed or warned in samples.

fails.and.warns <- qc %>%
                      select( sample, module, status ) %>%
                      filter( status %in% c( "WARN", "FAIL" ) ) %>%
                      arrange( sample )


##################################################################
### Summarizing Reports
##################################################################

# A summary and general statistics of the aggregated data.

# Summary of qc
summary( qc )

# Column names:
#
# module: fastqc modules
# nb_samples: the number of samples tested
# nb_pass, nb_fail, nb_warn: the number of samples that passed, failed and warned, respectively.
# failed, warned: the name of samples that failed and warned, respectively.


##################################################################
### General statistics
##################################################################

# The following table shows, for each sample, some general statistics such as the total number of reads, the length of reads, the percentage of GC content and the percentage of duplicate reads.


( fastQC.stats <- qc_stats( qc ) )

# Column names:
# pct.dup: the percentage of duplicate reads,
# pct.gc: the percentage of GC content,
# tot.seq: total sequences or the number of reads and
# seq.length: sequence length or the length of reads.


##################################################################
### Inspecting Problems
##################################################################

# Inspect problems to figure out what (if anything) is wrong with your data

# 1. R functions
# Used to inspect problems per either modules or samples using the following R functions

# Displays samples or modules that failed
qc_fails( qc )

# Displays samples or modules that warned
qc_warns( qc )

# Display which samples or modules that either failed or warned.
qc_problems( qc )

# 2. Input data

# Assess aggregated data from the previous table
qc

# 3. Output data

# Returns samples or FastQC modules with failures or warnings in a compact output format.
# For a stretched format, specify the argument compact = FALSE.

# The format and the interpretation of the outputs depend on the additional argument element, which value is one of c( “sample”, “module” ).

# If element = “sample” (default), results are samples with failed and/or warned modules.
# The results contain the following columns:

# sample (sample names);
# nb_problems (the number of modules with problems);
# module (the name of modules with problems).

#If element = “module”, results are modules that failed and/or warned in the most samples.
# The results contain the following columns:

# module (the name of module with problems);
# nb_problems (the number of samples with problems);
# sample (the name of samples with problems)

##################################################################
### Per Module Problems
##################################################################

# Modules that failed in the most samples:
qc_fails( qc, "module" )

# For each module, the number of problems (failures) and the name of samples, that failed, are shown below.

# Modules that warned in the most samples:
qc_warns( qc, "module" )

# For modules that either failed or warned:
qc_problems( qc, 'module' )

# The output above is in a compact format. For a stretched format, cancel the argument **compact**:

qc_problems( qc, "module", compact = FALSE )

# In the stretched format each row correspond to a unique sample and the status of each module is specified.
# To display problems for one or more specified modules change the argument **name**:
# Note that partial matching of name is allowed.
# For example, name = “Per sequence GC content” equates to name = “GC content”.

qc_problems( qc, "module",  name = "GC content" )


##################################################################
### Per Sample Problems
##################################################################

# Samples with one or more failed modules
qc_fails( qc, "sample" )

# For each sample, the number of problems (failures) and the name of modules, that failed, are shown.
# Samples with failed or warned modules:
# See which samples had one or more module with failure or warning
qc_problems( qc, "sample", compact = FALSE )

# To specify the name of a sample of interest, type this:
qc_problems( qc, "sample", name = "K562" )



# ##################################################################
# ### Building an HTML Report
# ##################################################################
# 
# # The function qc_report() can be used to build a report of FastQC outputs. It creates an HTML file containing FastQC reports of one or multiple samples.
# 
# # Inputs can be either a directory containing multiple FastQC reports or a single sample FastQC report.
# 
# ## Create a Multi-QC Report
# 
# #  Build a multi-qc report for the following demo QC directory:
# 
# # Build a report
# qc_report( qc.path = fastQC.dir,
#            result.file = wd.dir,
#            experiment = "Multi Quality Control for RNA-seq Datasets"
#            )
# 
# # An example of report is available at: <a href= “http://www.sthda.com/english/rpkgs/fastqcr/qc-reports/fastqcr-multi-qc-report.html”, target = "_blank"> fastqcr multi-qc report
# 
# ## Make a Report for each selected File
# 
# # Set path and file names
# 
# # Path to input files
# fastQC.files <- list.files( fastQC.dir, pattern = "*fastqc.zip$", full.names = TRUE)
# 
# # Path to the directory for output RNA-Seq single reports
# fastQC.sr.dir <- file.path( wd.dir, 'FastQC-SingleResults')
# fastQC.sr.files <- file.path( fastQC.sr.dir, substr( fastQC.files, 44, nchar( fastQC.files )-10 ) )
# fastQC.sr.name <- substr( fastQC.files, 44, nchar( fastQC.files )-11 )
# 
# # One sample QC report with interpretation
# 
# ##################################################################
# ### Importing and Plotting a FastqQC QC Report
# ##################################################################
# 
# for ( i in 1:length( fastQC.files ) ) {
#   qc_report( qc.path = fastQC.files[ i ],
#              result.file = 'Sample_Report',
#              interpret = TRUE )
# }
# 
# 
# 
# # Elements contained in the qc object
# names( qc[[1]] )
# 
# # The function qc_plot() is used to visualized the data of a specified module.
# # Allowed values for the argument modules include one or the combination of:
# 
# # “Summary”,
# # “Basic Statistics”,
# # “Per base sequence quality”,
# # “Per sequence quality scores”,
# # “Per base sequence content”,
# # “Per sequence GC content”,
# # “Per base N content”,
# # “Sequence Length Distribution”,
# # “Sequence Duplication Levels”,
# # “Overrepresented sequences”,
# # “Adapter Content”
# 
# par( mfrow = c( 2, 2 ) )
# qc_plot( qc, "Per base sequence quality" )
# qc_plot( qc, "Per sequence quality scores" )
# qc_plot( qc, "Per base sequence content" )
# qc_plot( qc, "sequence length distribution" )
# 





##################################################################
### Software and Code Used
##################################################################
sessionInfo()
