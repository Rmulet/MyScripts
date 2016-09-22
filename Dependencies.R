#!/usr/bin/Rscript

# CRAN REPOSITORY
install.packages("PopGenome",repos="http://cran.r-project.org")
install.packages("stringr",repos="http://cran.r-project.org")

# BIOCONDUCTOR:

source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite(TxDb.Hsapiens.UCSC.hg19.knownGene)