#!/usr/bin/Rscript

# CRAN REPOSITORY
if (length(setdiff(c("PopGenome","stringr"), rownames(installed.packages()))) > 0) {
  install.packages(setdiff(c("PopGenome","stringr"), rownames(installed.packages()),repos="http://cran.r-project.org"),INSTALL_opts = c('--no-lock'))
}

# BIOCONDUCTOR:

if (length(setdiff(c("GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene","Biostrings"), rownames(installed.packages()))) > 0) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(setdiff(c("GenomicRanges","TxDb.Hsapiens.UCSC.hg19.knownGene","Biostrings"), rownames(installed.packages())))
}
  
# PREPARE GENE LIST

setwd("~/Documents/2_GenomicsData/Final/GeneByGene")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels0(txdb) # Resets seqlevels
seqlevels(txdb) <- sapply(1:22,function(x){paste("chr",x,sep="",collapse="")})  
grgenes <- genes(txdb)
genestable <- data.frame(name=grgenes$"gene_id",chr=seqnames(grgenes),start=start(grgenes),end=end(grgenes))
save(genestable,file="GenesTable.RData")