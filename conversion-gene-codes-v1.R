# Institution: The University of Melbourne
# Author: Marc Rigau
# Date: 24/1/2022
# Last update: 21/4/2022
# Version: 1
#################################################################

# Source: https://bioconductor.riken.jp/packages/3.7/bioc/vignettes/biomaRt/inst/doc/biomaRt.html

# Notes

# The biomaRt package, provides an interface to a growing collection of databases implementing the BioMart software suite. The package enables retrieval of large amounts of data in a uniform way without the need to know the underlying database schemas or write complex SQL queries.
# Examples of BioMart databases are Ensembl, Uniprot and HapMap.
# These major databases give biomaRt users direct access to a diverse set of data and enable a wide range of powerful online queries from R.

#################################################################

# Load specific library from Bioconductor
# BiocManager::install(c("biomaRt"))
library("biomaRt")

# Set a proxy if runing into trouble with useMart()
#Sys.setenv("http_proxy" = "http://my.proxy.org:9999")

#################################################################
# Selecting a BioMart database and dataset

listMarts()
useMart( "ensembl" )

# Connect to a specific BioMart database
ensembl=useMart("ensembl")

# View all datasets in the database
listDatasets(ensembl)

# If known dataset, select a BioMart database in a single step:
ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)

# Alternatively
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")


#################################################################
# Build a biomaRt query

# Prepare the three main parameters for a BiomaRT querry

# 1 Filters
filters = listFilters(ensembl)
filters[1:5,]

# Attributes
attributes = listAttributes(ensembl)
attributes[1:5,]

# Main BiomaRT command:
# attributes: is a vector of attributes that one wants to retrieve (= the output of the query).
# filters: is a vector of filters that one wil use as input to the query.
# values: a vector of values for the filters. In case multple filters are in use, the values argument requires a list of values where each position in the list corresponds to the position of the filters in the filters argument (see examples below).
# mart: is an object of class Mart, which is created by the useMart() function.

# Note: for some frequently used queries to Ensembl, wrapper functions are available: getGene() and getSequence().
# These functions call the getBM() function with hard coded filter and attribute names.

# BioMart query

ensembl_ids_samples = c("202763_at","209310_s_at","207500_at")

getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
      filters = 'ensembl_gene_id', 
      values = affyids, 
      mart = ensembl)






















