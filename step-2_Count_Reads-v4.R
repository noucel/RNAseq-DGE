##################################################################
### Libraries
##################################################################
library( tidyverse ) # Kernel for statistics in science fields.
library( Rsubread ) # Map the reads to the reference genome (eg. the human genome “hg19”) (Y et al. 2013).

##################################################################
### Count Reads for each Feature
##################################################################
# Assign mapped sequencing reads to genes based on the index previously built.

# Set the working space directory
ws.dir <- "/Volumes/WolfeCreek/Marc/RNAseq"

# Define the directories for the alignment files
mapping.dir <- paste0( ws.dir, "/Mapping" )
ref.genome.dir <- "~/UniMelb/Datasets/RNAseq/GENCODE-Human_Release-40-Index"

counts <- list.files( mapping.dir, pattern = ".bam$" )
counts

# Assign mapped sequencing reads to specified genomic features
fc <- featureCounts( file.path( mapping.dir, counts ),
                     annot.ext = file.path( ref.genome.dir, "gencode_v40_annotation.gtf" ),
                     isGTFAnnotationFile = TRUE,
                     isPairedEnd = TRUE,
                     GTF.attrType = "gene_name",
                     nthreads = 64,
                     requireBothEndsMapped = TRUE,
                     verbose = FALSE
)

# Summary of the output object
summary( fc )
dim( fc$counts )

# Name by unique and stable GEO accession number
colnames( fc$counts ) <- substr( colnames( fc$counts ), 13, 22 )
fc$counts[ 1:2, ]

# Define the directory for the count files
ws.dir <- "~/UniMelb/Datasets/RNAseq"
counts.dir <- paste0( ws.dir, "/Counts" )

# Save the output data in a tabular separated table format
write.table(fc$counts,
            file = file.path( counts.dir, "featureCounts_Homo_sapiens_RNA-Seq.txt" ),
            sep="\t", quote=F, append=F, row.names = T
            )


##################################################################
### Quality Control and Stats
##################################################################

# Import or built an experiment design table to assign sample information.

GEO <- substr( counts, 13, 22 ) # Name by unique and stable GEO accession number
Species <- c( rep( 'Homo Sapiens', length( counts ) ) )
Speciment <- c( rep( 'U266', 6), rep( 'K562', 3), rep( 'U266', 2), rep( 'K562', 2) )
Batch <- c( '4h', '24h', '4h', '24h', '4h', '24h', 'Cnt', 'Cnt', 'Cnt', 'scrbl', 'scrbl', 'WT', 'WT'  )
Replicate <- c( 1, 1, 2, 2, 3, 3, 1, 2, 3, 1, 2, 1, 2 )

# Assemble a design dataframe
sample.annotation <- as.data.frame(  cbind(Species, Speciment, Batch, Replicate), row.names = GEO )

# List samples in an object
samples <- as.character( rownames( sample.annotation) )

# Grup speciments for the analysis question
group <- factor( sample.annotation$Speciment )
table( group )

# Comon name samples
tags <- as.factor( paste( group, Batch, Replicate ) )
tags


##################################################################
### Basic Quality Control Plots
##################################################################

# Define the output for Rplots directory
Rplots.dir <- paste0( ws.dir, "/Rplots" )

# Density plots of log-intensity distribution for each library.
x.max <- max( log( fc$counts, 10 ) )

# Density plot of raw read counts (log10)
png( file.path( Rplots.dir, 'raw_read_counts_per_gene.density.png' ) )
for ( i in 1:length( samples ) ) {
  logcounts <- log( fc$counts[ , i ], 10 )
  d <- density( logcounts )
  if ( i == 1 ) plot( d, xlim = c( 0, x.max), ylim = c( 0, 0.24 ), main = "", xlab = "Raw read counts per gene (log10)", ylab = "Density", col = colors()[ 10 ], lwd = 3 )
  else lines( d, col = colors()[ 5 * i ], lwd = 3 )
}
dev.off()


# Boxplots of the raw read counts after log10 transformation
png( file.path( Rplots.dir, 'raw_read_counts_per_gene.boxplot.png' ) )
logcounts <- log( fc$counts, 10 )
# Ignore -Inf values; instead consider non-expressed gene 'NA' or replace it by a constant EPSILON
EPSILON = NA # 1E-10
logcounts[ is.infinite( logcounts ) ] <- EPSILON
colnames( logcounts ) <- tags
boxplot( logcounts, main = "", alpha= 0.1, xlab = "", ylab = "Raw read counts per gene (log10)" )
#boxplot( logcounts, col = colors()[ seq( 10, 10*nlevels( common.names ), by = 10 ) ], main = "", xlab = "", ylab = "Raw read counts per gene (log10)" )
dev.off()


# Cluster hierarchicaly to investigate the relationship between samples.
# Make a heatmap of the highest expressed genes to euclide distances.

# Select data for the most highly expressed genes
threshold <- 100
temp = order( rowMeans( fc$counts ), decreasing = TRUE )[ 1:threshold ]
top.gene.counts <- fc$counts[ temp, ]

# To better understand what biological effect lies under this clustering, use the samples annotation for labeling.
colnames( top.gene.counts ) <- tags

# Heatmap
png( file.path( Rplots.dir, 'high_expr_genes.heatmap.png' ) )
heatmap( top.gene.counts, col = heat.colors( threshold, alpha= 1, rev = FALSE ) )
dev.off()



##################################################################
### Principal Component Analysis
##################################################################

# Principal Component Analysis (PCA) from the cmdscale function (from the stats package) performs a classical multidimensional scaling of a data matrix.
# Reads counts need to be transposed before being analysed with the cmdscale functions.
# i.e. genes should be in columns and samples should be in rows.

# Select data for the most highly expressed genes
threshold <- 100
highly.expressed.gene.counts = order( rowMeans( fc$counts ), decreasing = TRUE )[ 1:threshold ]
highly.expressed.gene.counts <- fc$counts[ highly.expressed.gene.counts, ]
# Annotate the data with condition group as labels
colnames( highly.expressed.gene.counts ) <- tags
# Transpose the data to have variables (genes) as columns
data.for.PCA <- t( highly.expressed.gene.counts )
dim( data.for.PCA )

# The cmdscale function calculates a matrix of dissimilarities from your transposed data.
# Information about the proportion of explained variance is also provided by calculating Eigen values.

# Calculate matrix of dissimilarities (MDS)
# Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis (Gower, 1966).
mds <- cmdscale( dist( data.for.PCA ), k = 3, eig = TRUE )  
# k = the maximum dimension of the space which the data are to be represented in
# eig = indicates whether eigenvalues should be returned

# The variable mds$eig provides the Eigen values for the first 8 principal components:
mds$eig

# Plotting this variable as a percentage determines how many components can explain the variability in the dataset and how many dimensions you should be looking at.
# Transform the Eigen values into percentage
eig.percentage <- mds$eig * 100 / sum( mds$eig )

# plot the PCA
png( file.path( Rplots.dir, 'PCA_PropExplainedVariance.png' ) )
barplot(eig.percentage,
        las = 1,
        xlab = "Dimensions", 
        ylab = "Proportion of explained variance (%)", y.axis = NULL,
        col = "skyblue" )
dev.off()

# In most cases, the first 2 components explain more than half the variability in the dataset and can be used for plotting.
# The cmdscale function run with default parameters performs a PCA on the given data matrix.
# The plot function provides scatter plots for individuals representation.

# Calculate MDS
mds <- cmdscale( dist( data.for.PCA ) ) # Performs MDS analysis 

# Samples representation
png( file.path( Rplots.dir, 'PCA_Dim1vsDim2.png' ) )
plot( mds[ , 1 ], -mds[ , 2 ], type = "n", xlab = "Dimension 1", ylab = "Dimension 2", main = "")
text( mds[ , 1 ], -mds[ , 2 ], rownames( mds ), cex = 0.8) 
dev.off()

# The PCA plot of the first two components should show a clear separation of sample groups across the 1st or 2nd dimension.
# Each group of samples represent a cluster.
# Different factors can show other type of variability among the data, often dirven by biological biais such as treatment, disease, or cell type.



##################################################################
### Software and Code Used
##################################################################
sessionInfo()



