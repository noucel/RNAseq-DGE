##################################################################
### Differential Gene Expression Analysis
##################################################################
# Load packages
library( tidyverse )
library( DT )
library( edgeR )
library( Glimma )

##################################################################

# Read Count directory
read.counts.dir <- "~/UniMelb/Datasets/RNAseq/Counts"

feature.counts <- file.path( read.counts.dir, 'featureCounts_Homo_sapiens_RNA-Seq.txt' )

mycounts <- read.table( feature.counts )

# The following command is to create a DGE list, which is not part of this workflow (other guide)
# Create a differential gene expression (DGE) list
y <- DGEList( fc, lib.size = colSums( fc$counts ),
              norm.factors = calcNormFactors( fc$counts ),
              samples = sample.common.names,
              group = group
              )




# Make sample annotation table

Speciment <- c( rep( 'U266', 2), rep( 'K562', 2) )
Replicate <- c( 1, 2, 1, 2 )
GEO <- substr( sample.name, 13, 22 ) # Name by unique and stable GEO accession number


sample.anno <- as.data.frame( cbind( Speciment, Replicate, GEO ) )
sample.anno

# Create a dataframe with annotations for the differential gene expression analysis (DGEList).
dge <- readDGE( files = 'featureCounts_Homo_sapiens_RNA-Seq.txt read.counts.dir',
                path = read.counts.dir
                group = sample.anno$Speciment, labels=sample.anno$GEO, header = F )

m <- as.matrix.DGEList( fc  )





rm(dge)



# Obtain the expression counts from previous alignment step
counts <- fc$counts
dim( counts )
head( counts )

# All counts in a glance
hist(
  log( rowMeans( counts ) ),
  breaks = 100,
  main = "Count Distribution",
  xlab= 'Count per sequence (ln)'
)





# The following command is to create a DGE list, which is not part of this workflow (other guide)
# # Create a differential gene expression (DGE) list
# y <- DGEList( fc$counts, lib.size = colSums( fc$counts ),
#               norm.factors = calcNormFactors( fc$counts ),
#               samples = samples$samplename,
#               group = samples$condition
#               )


# Filtering
# Keep genes with least 'x' count-per-million reads (cpm) in at least 'n' samples
isexpr <- rowSums( cpm( gene.counts ) > 1) >= 4
table( isexpr )

gene.counts <- gene.counts[ isexpr, ]
genes <- rownames( gene.counts )

dim( gene.counts )

# The limma package (since version 3.16.0) offers the voom function that normalises read counts and applies a linear model to the normalised data before computing moderated t-statistics of differential expression.
# limma was originally developed for differential expression analysis of microarray data.
# Load required libraries
#library( limma )

# Check samples grouping
experiment.design.ord[ colnames( gene.counts ), ]$Speciment == group

# Create design matrix for limma
design <- model.matrix( ~0 + group )

# Substitute "group" from the design column names
colnames( design ) <- gsub( "group", "", colnames( design ) )

# Check your design matrix
design

# Calculate normalization factors between libraries
nf <- calcNormFactors( gene.counts )

# Normalise the read counts with 'voom' function
# voom is a function in the limma package that modifies RNA-Seq data for use with limma.
y <- voom( gene.counts, design, lib.size = colSums( gene.counts ) * nf )


# Extract the normalised read counts
counts.voom <- y$E

# Save normalised expression data into output dir
write.table( counts.voom,
             file = file.path( read.counts.dir, "featureCounts_Homo_sapiens_RNA-Seq.txt" ),
             row.names = T , quote = F, sep = "\t" )

# Fit linear model for each gene given a series of libraries
fit <- lmFit( y, design )

# Construct the contrast matrix corresponding to specified contrasts of a set of parameters
cont.matrix <- makeContrasts( U266-K562, levels = design )
cont.matrix 

# Compute estimated coefficients and standard errors for a given set of contrasts
fit <- contrasts.fit( fit, cont.matrix )

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit <- eBayes( fit )
options( digits = 3 )

# check the output fit
dim( fit )

# The topTable function summarises the output from limma in a table format.
# Significant differential expression (DE) genes for a particular comparison can be identified by selecting genes with a p-value smaller than a chosen cut-off value and/or a fold change greater than a chosen value in this table.
# By default the table will be sorted by increasing adjusted p-value, showing the most significant DE genes at the top.

# Set adjusted pvalue threshold and log fold change threshold
mypval = 0.01
myfc = 3

# Obtain the coefficient name for the comparison  of interest
colnames( fit$coefficients )

# Obtain the output table for the 10 most significant DE genes for this comparison
mycoef = colnames( fit$coefficients )
topTable( fit, coef = mycoef )

# Obtain the full table ("n = number of genes in the fit")
limma.res <- topTable( fit, coef = mycoef, n = dim( fit )[ 1 ] )

# Obtain significant DE genes only (adjusted p-value < mypval)
limma.res.pval <- topTable( fit, coef = mycoef, n = dim( fit )[ 1 ], p.val = mypval )
dim( limma.res.pval )

# Obtain significant DE genes with low adjusted p-value high fold change
limma.res.pval.FC <- limma.res.pval[ which( abs( limma.res.pval$logFC ) > myfc ), ]
dim( limma.res.pval.FC )

# Write limma output table for significant genes into a tab delimited file
filename.dge = paste( "DGE_", mycoef, "_pvalue-", mypval, "_logFoldChange-", myfc, ".txt", sep="" )
filename.dge.path = file.path( read.counts.dir, filename.dge )

write.table( limma.res.pval.FC, file=filename.dge.path, row.names = T, quote = F, sep = "\t" )


# Plot a specified coefficient of the linear model
volcanoplot(fit = fit,
            coef = 1,
            style = "p-value",
            highlight = 100,
            names = rownames( fit ),
            hl.col = "blue",
            xlab = "Log2 Fold Change"
)


# Oorganize the labels nicely using the "ggrepel" package and the geom_text_repel() function
library(ggrepel)

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they have a Higher- or Lower- expression (log2FoldChange respectively positive or negative)

# add a column of NAs
limma.res.pval.FC$diffexpressed <- "Similar"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
limma.res.pval.FC$diffexpressed[limma.res.pval.FC$logFC > 0.6 & limma.res.pval.FC$P.Value < 0.05] <- "Higher"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
limma.res.pval.FC$diffexpressed[limma.res.pval.FC$logFC < -0.6 & limma.res.pval.FC$P.Value < 0.05] <- "Lower"


# ggplot( data=limma.res.pval.FC, aes( x = limma.res.pval.FC$logFC, y = -log10( limma.res.pval.FC$P.Value ), col = limma.res.pval.FC$diffexpressed, label = rownames( limma.res.pval.FC ) ) ) +
#   geom_point() + 
#   theme_minimal() +
#   geom_text_repel() +
#   scale_color_manual(values=c("blue", "black", "red")) +
#   geom_vline(xintercept=c(-0.6, 0.6), col="red") +
#   geom_hline(yintercept=-log10(0.05), col="red")




##################################################################
### Gene Annotation
##################################################################

# Add annotations of the BioMart database to obtain information about the target genes.
# Load the Ensembl annotation for human genome form the BioConductor repository
# BiocManager::install("biomaRt")
library( biomaRt )

# Select a BioMart database and dataset
listMarts()
useMart( "ensembl" )

# Connect to a specific BioMart database
ensembl = useMart( "ensembl" )

# View all datasets in the database
listDatasets( ensembl )

# If known dataset, select a BioMart database in a single step:
ensembl = useDataset( "hsapiens_gene_ensembl", mart = ensembl )

# Alternatively, if the database is already known use the following command
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")


#################################################################
# Build a biomaRt query

# Prepare the three main parameters for a BiomaRT querry

# 1 Filters
filters = listFilters( ensembl )
filters[ 1:10,  ]

# Attributes
attributes = listAttributes( ensembl )
attributes[ 1:25, ]

# Main BiomaRT command:
# attributes: is a vector of attributes that one wants to retrieve (= the output of the query).
# filters: is a vector of filters that one wil use as input to the query.
# values: a vector of values for the filters. In case multple filters are in use, the values argument requires a list of values where each position in the list corresponds to the position of the filters in the filters argument (see examples below).
# mart: is an object of class Mart, which is created by the useMart() function.

# Note: for some frequently used queries to Ensembl, wrapper functions are available: getGene() and getSequence().
# These functions call the getBM() function with hard coded filter and attribute names.

# BioMart query

# Pull gene IDs from limma output table
hgnc_genes <- as.character( rownames( limma.res.pval.FC ) )
length( hgnc_genes )

detags.IDs <- getBM( attributes = c('hgnc_symbol', 'entrezgene_id', 'ensembl_gene_id', 'ensembl_transcript_id', 'chromosome_name', 'band', 'strand', 'description' ),
                     filters = 'hgnc_symbol', 
                     values = hgnc_genes, 
                     mart = ensembl )
dim( detags.IDs )
head( detags.IDs )


# Remove duplicates
detags.IDs.matrix <- detags.IDs[ -which( duplicated( detags.IDs$hgnc_symbol ) ), ]

# Select genes of interest only
rownames( detags.IDs.matrix ) <- detags.IDs.matrix$hgnc_symbol
hgnc_genes.annot <- detags.IDs.matrix[ as.character( hgnc_genes ), ]

# Join the two tables
limma.res.pval.FC.annot <- cbind( hgnc_genes.annot, limma.res.pval.FC )

# Check the annotated table
head( limma.res.pval.FC.annot )

# Write limma output table for significant genes into a tab delimited file
filename.dge.annot = paste( substr( filename.dge, 1, nchar( filename.dge) -4 ), '_Annotation', ".txt", sep="" )
filename.dge.annot.path = file.path( read.counts.dir, filename.dge.annot )

write.table( limma.res.pval.FC.annot, file=filename.dge.annot.path, row.names = T, quote = F, sep = "\t" )















##################################################################
### Software and Code Used
##################################################################
sessionInfo()
