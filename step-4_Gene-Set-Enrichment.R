##################################################################
### Gene Set Enrichment
##################################################################

# Gene Ontology (GO) enrichment is a method for investigating sets of genes using the Gene Ontology system of classification, in which genes are assigned to a particular set of terms for three major domains: cellular component, biological process and molecular function.
# The GOstats package can test for both over and under representation of GO terms using the Hypergeometric test.
# The output of the analysis is typically a ranked list of GO terms, each associated with a p-value.
# The Hypergeometric test will require both a list of selected genes (i.e. your DE genes) and a “universe” list (e.g. all genes annotated in the genome you are working with), all represented by their “EntrezGene” ID.

# load the library
# BiocManager::install( "GOstats" )
# BiocManager::install( "org.Hs.eg.db" )
library( GOstats )
library( org.Hs.eg.db )




# fc.entrezgene <- featureCounts( file.path( mapping.dir, count.files ),
#                                 annot.ext = file.path( reference.genome.dir, "gencode_v40_annotation.gtf" ),
#                                 isGTFAnnotationFile = TRUE,
#                                 isPairedEnd = TRUE,
#                                 GTF.attrType = "gene_id", # Gene IDs are Ensembl's gene ID begging with ENS for Ensembl, and then a G for gene. 
#                                 nthreads = 8,
#                                 requireBothEndsMapped = TRUE,
#                                 verbose = FALSE
# )
# 
# # Name by unique and stable GEO accession number
# colnames( fc.entrezgene$counts ) <- substr( colnames( fc.entrezgene$counts ), 13, 22 )
# #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

top.genes.annot <- read.table( filename.dge.annot.path, row.names = T )

# Define list of genes of interest (DE genes - HUGO Gene Nomenclature Committee)
entrezgeneids <- as.character( top.genes.annot$ensembl_gene_id )
length( entrezgeneids )

# Define the universe
universeids <- rownames( mycounts )
length( universeids )


# Before running the hypergeometric test with the hyperGTest function, the parameters for the test (gene lists, ontology, test direction) have to  be defined, as well as the annotation database to be used.
# The ontology to be tested can be any of the three GO domains: biological process (“BP”), cellular component (“CC”) or molecular function (“MF”).

# Define the p-value cut off for the hypergeometric test
hgCutoff <- 0.05
params <- new( "GOHyperGParams",
               annotation = "org.Hs.eg",
               geneIds = entrezgeneids,
               universeGeneIds = universeids,
               ontology = "BP",
               pvalueCutoff = hgCutoff,
               testDirection = "over"
)

#  Run the test
hg <- hyperGTest( params )

# Check results
hg