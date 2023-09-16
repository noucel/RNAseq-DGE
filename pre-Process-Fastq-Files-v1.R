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

# https://anilchalisey.github.io/rseqR/

# Install dependences:
# Fastqc (installed in the previuos quality control file)
# MultiQC



# Note: Once all dependencies are installed, then rseqR may be installed as follows:
# Install rseqR
devtools::install_github("anilchalisey/rseqr", build_vignettes = TRUE)

# Trim fastq files

##################################################################
#### Pre-process the Raw Sequencing Data
##################################################################






##################################################################
### Software and Code Used
##################################################################
sessionInfo()
