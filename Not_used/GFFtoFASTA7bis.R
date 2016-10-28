#!/usr/bin/Rscript

# setwd("~/Documents/2_GenomicsData/Final/GeneByGene")
# gffpath <- "chr21.gff"; ref.chr <- "chr21.fa"; chrom <- "21"

# This script creates a FASTA-like file with numeric codes identifying each feature. Thus, 9 is a gene,
# 8 is a CDS, 5 is 5'-UTR and so on. For most features, it takes them directly from the GFF; for the coding
# regions, it applies check_fold to determine their degeneracy (4 for 4-fold, 0 for 0-fold)

# NOTE: The GFF file is 1-based, like R, i.e. one base is denoted as 1:1.
# WARNING: The v8 version is intended for use of UCSC GFF files without gene notation.


suppressMessages(library(PopGenome))
suppressMessages(library(stringr))

#############################
## IMPORT AND PREPARE DATA ##
#############################

args <- commandArgs(trailingOnly = TRUE) # Import arguments from command line
gffpath <- args[1] # Name of GFF file of the chromosome
ref.chr <- args[2] # Name of the reference chromosome in FASTA

chrom <- str_match(ref.chr,"chr([0-9]+).")[2]

if (args[1] == "-h" | args[1] == "--help") {
  cat("\nGFFtoFASTA - A script that creates a FASTA-like file with numeric codes identifying each feature.\n")
  cat("\nUsage: GFFtoFASTA.R [CHR GFF FILE] [CHR FASTA FILE]\n\n")
  quit()
}

#########################
## FUNCTION DEFINITION ##
#########################

# PopGenome:::codonise64() to call codonise. Remember it takes a one-row matrix as argument.

T <- codontable()
trips <- T$Triplets
code <- T$Protein[1,]   

# codontable() # 1st matrix - aa of the triplets, each row is a different code (first standard)
# 2nd matrix - nucleotide combinations of each aa (row), matching the matrix above. U=1, C=2, A=3, G=4.
# region@region.data@CodingSNPS - Coding SNPs TRUE/FALSE

check_fold <- function(triplet,pos) {
  # CHECK THAT NO ARGUMENTS ARE MISSING
  if (missing(triplet) | missing(pos)){
    stop("Arguments are missing")
  }
  # EVALUATE INPUT VARIABLES:
  comb <- which(apply(trips,1,function(x){identical(x[1:3],triplet)}))
  if (comb == 65) {return(fold=NA)} # For the 555 case
  aa <- code[comb] # Original aminoacid of the triplet
  # EVALUATE DEGENERACY:
  fold <- 0
  for (i in 1:4) {
    triplet[pos+1] <- i # Pos in cod.pos is 0/1/2, but indices in R start at 1
    comb <- which(apply(trips,1,function(x){identical(x[1:3],triplet)}))
    altaa <- code[comb]
    # cat(c(triplet,",",altaa,"\n"))
    if (altaa == aa) {fold <- fold+1}
  }
  return (fold)
}

# Associative array: fold for each possible position
assocmat <- t(apply(trips,1,function(x){sapply(0:2,function(y){check_fold(x,y)})}))
assocmat <- assocmat[-nrow(assocmat),]
assocmat[assocmat == 1] <- 0

#####################################
## ASSIGNATION OF SEQUENCE CLASSES ##
#####################################

## RETRIEVING SEQUENCE CLASSES TO VECTOR

suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
library(GenomicRanges)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb) <- sprintf("chr%s",chrom)
genes <- transcripts(txdb) # The function genes() returns HUGO genes, not UCSC transcripts
CDS <- cds(txdb)
exons <- exons(txdb)
fiveUTR <- unlist(fiveUTRsByTranscript(txdb))
threeUTR <- unlist(threeUTRsByTranscript(txdb))

file.info <- .Call("get_dim_fasta",ref.chr) # .Call invokes the C function get_ind_fasta (src folder): start and end positions of chr.
gffseq <- as.integer(rep(NA,file.info[[1]][2])) # Stores the Integer is smaller than numeric
# We assign features/classes from less to more conserved: 0 > 5 > 3 > Intron (9) > 4 > Intergenic (NA).
# Non-assigned genes (9) will be introns. Exons can be coding (CDS) or non-coding (8).
# CDS (7) will be replaced with 0/4 later.
features <- c("genes"=9L,"exons"=8L,"CDS"=7L,"threeUTR"=3L,"fiveUTR"=5L)
for (n in 1:length(features)) {
  feature <- eval(as.symbol(names(features[n])))
  symbol <- features[n]
  for (i in 1:length(feature)) {
    gffseq[start(feature[i]):end(feature[i])] <- symbol
  }
}
sum(gffseq == 9,na.rm=T) # Sanity check. 18950868 gene.

# Alternatively, we could use the 'TxDb.Hsapiens.UCSC.hg19.knownGene' package containing all regions in GRanges objects.
# However, importing from the GFF allows for greater control and it is also closer to what PopGenome uses.

####################################
## DEGENERACY IN CODING SEQUENCES ##
####################################

## LOAD GFF FILE

gff.table   <- read.table(gffpath,sep="\t",colClasses=c(rep("character",3),
                                                        rep("integer",2),rep("character",2),"character","character")) # Replace with chr!!!
cdsgff <- gff.table[gff.table[,3] == "CDS",]
cdsgff[,9] <- str_match(cdsgff[,9],"gene_id ([:alnum:]+)")[,2] # Keep the names
cdsgff <- cdsgff[,c(1,3,4,5,7,8,9)]
colnames(gff.table) <- c("chr","feature","start","end","strand","frame","ID")

## OBTAIN CODING SEQUENCE AND CHECK FOLD 

CHR <- .Call("get_ind_fasta",ref.chr,1,file.info[[1]][2]) # .Call invokes get_ind_fasta, which extracts the letters in FASTA to numbers.

komplement <- c(4,3,2,1,5) # For regions in the reverse strand: opposite of T=1,C=2,G=3,A=4 (vectorisation)
for (x in unique(cdsgff$ID)) { # x for each gene/transcript in the GFF file
  tab <- cdsgff[cdsgff$ID == x,][,3:5]
  if(tab$strand[1] == "-") { tab <- tab[,c(2:1,3)]} # For negative strand, gene goes from END to START
  # Concatenate the CDS corresponding to each gene:
  trans <- as.vector(unlist(apply(tab,1,function(x){CHR[x[1]:x[2]]}))) # Translation vector
  if(tab$strand[1] == "-") { trans <- komplement[trans] }

  # Extract the codons/triplets of the concatenated coding region:
  fold <-vector() # Fold degeneracy of the entire gene
  a <- matrix(trans,ncol=3,byrow=TRUE)
  indices <- seq(0,length(trans)-1)%%3+1 # +1 because we use assocmat, -1 because we start at 0
  for (i in 1:length(trans)){ # Find degeneracy in every position
    mpos <- ceiling(i/3) # Horizontal position in the matrix (i.e. triplet number)
    ncod <- as.numeric(PopGenome:::codonise64(a[mpos,,drop=F]))
    # test <- c(test,code[ncod])
    fold <- c(fold,assocmat[ncod,indices[i]]) # Add the position corresponding to that region
  }
  # Note that, for genes in the negative strand, the coordinates go from big to small (i.e. 5' to 3' in that strand).
  # Thus, since the coordinates are already switched, the 'fold' vector will be copied in the proper order.
  
  # Fraction the 'fold' vector according to the CDS coordinates.
  pos1 <- 1  
  for (n in 1:nrow(tab)){
    pos2 <- pos1+abs(tab[n,2]-tab[n,1]) # Position in the fold object (concatenated CDS). Add +1 because it starts at pos1!
    local.fold <- fold[pos1:pos2]
    gffseq[tab[n,1]:tab[n,2]][local.fold == 0] <- 0
    # Four-fold are less constrained that 5UTR, 3UTR and introns:
    relaxed <- gffseq[tab[n,1]:tab[n,2]] == 8 | gffseq[tab[n,1]:tab[n,2]] == 7 # Relaxed positions (i.e. not 5, 3, 0 or 9)
    gffseq[tab[n,1]:tab[n,2]][relaxed & local.fold == 4] <- 4 # Hence, they can be replaced with 4-fold.
    pos1 <- pos2+1 
  }
  print(x) # ID of the analysed gene
}
Sys.time()-init

save(gffseq,file=sprintf("gffseq_chr%s.RData",chrom))
