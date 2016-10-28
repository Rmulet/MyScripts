# setwd("~/Documents/2_GenomicsData/TestPopGenome")
# ref.chr = "chr22.gff"

# This script creates a FASTA-like file with numeric codes identifying each feature. Thus, 9 is a gene,
# 8 is a CDS, 5 is 5'-UTR and so on. For most features, it takes them directly from the GFF; for the coding
# regions, it applies check_fold to determine their degeneracy (4 for 4-fold, 0 for 0-fold)

library(PopGenome)
library(stringr)

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

## VECTOR CONTAINING THE GFF CLASSES
file.info <- .Call("get_dim_fasta",ref.chr) # .Call invokes the C function get_ind_fasta (src folder): start and end positions of chr.
gffseq1 <- as.integer(rep(NA,file.info[[1]][2])) # Stores the Integer is smaller than numeric

# Before applying the loops that do not depend on strand, reduce complexity with ranges(genes)

## RETRIEVING SEQUENCE CLASSES

init <- Sys.time()
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb,force=TRUE) <- "chr22" # Restricts to chr22.
init <- Sys.time()
genes <- transcripts(txdb) # The function genes() returns HUGO genes, not UCSC transcripts
cds <- cds(txdb)
fiveUTR <- unlist(fiveUTRsByTranscript(txdb))
threeUTR <- unlist(threeUTRsByTranscript(txdb))
features <- c("genes"=9L,"cds"=8L,"fiveUTR"=5L,"threeUTR"=3L)
# init <- Sys.time()
for (n in 1:length(features)) {
  feature <- eval(as.symbol(names(features[n])))
  symbol <- features[n]
  for (i in 1:length(feature)) {
    gffseq1[start(feature[i]):end(feature[i])] <- symbol
  }
}
Sys.time()-init

sum(gffseq == 9,na.rm=T) # 19292012

## LOAD GFF FILE

init <- Sys.time()
gff.table   <- read.table("chr22.gff",sep="\t",colClasses=c(rep("character",3),
                rep("integer",2),rep("character",2),"character","character")) # Replace with chr!!!
init <- Sys.time()
gff.table[,9] <- str_match(gff.table[,9],"ID=(?:[\\w_]*:)?([a-z-0-9.]+)")[,2] # Keep the names
Sys.time()-init
gff.table <- gff.table[,c(1,3,4,5,7,8,9)]
colnames(gff.table) <- c("chr","feature","start","end","strand","frame","ID")
cdsgff <- gff.table[gff.table$feature == "CDS",]
Sys.time()-init

gffseq <- as.integer(rep(NA,file.info[[1]][2])) # Stores the Integer is smaller than numeric
# We assign features from less to more conserved: 0 > 5 > 3 > Intron (9) > 4 > Intergenic (NA).
# Non-assigned genes (9) will be introns. Exons can be coding (CDS) or non-coding (8).
# CDS (7) will be replaced with 0/4 later.
init <- Sys.time()
features <- c("gene"=9L,exon="8L","CDS"=7L,"three_prime_UTR"=3L,"five_prime_UTR"=5L)
for (n in 1:length(features)) {
  feature <- gff.table[gff.table$feature == names(features[n]),]
  symbol <- features[n]
  for (i in 1:nrow(feature)) {
    gffseq[feature[i,3]:feature[i,4]] <- symbol # 19353446
  }
}
sum(gffseq == 9,na.rm=T) # 19353446

Sys.time()-init

####################################
## DEGENERACY IN CODING SEQUENCES ##
####################################

## OBTAIN CODING SEQUENCE AND CHECK FOLD 

CHR <- .Call("get_ind_fasta",ref.chr,1,file.info[[1]][2]) # .Call invokes get_ind_fasta, which extracts the letters in FASTA to numbers.
strand(cds) <-"-" # Select only a given strand

init <- Sys.time()
# sapply(unique(cdsgff$ID),function(x) {
for (x in unique(cdsgff$ID)) {
  tab <- cdsgff[cdsgff$ID == x,][,3:5]
  fold <-vector()
  if(tab$strand[1] == "+") {
    trans  <- as.vector(unlist(apply(tab,1,function(x){CHR[x[1]:x[2]]})))
  } else if (tab$strand[1] == "-") { # It is ordered from start (higher) to end (lower)
    trans <- as.vector(unlist(apply(tab,1,function(x){CHR[x[2]:x[1]]})))
    trans <- komplement[trans]
  }
  a <- matrix(trans,ncol=3,byrow=TRUE)
  indices <- seq(0,length(trans)-1)%%3+1 # +1 because we use assocmat, -1 because we start at 0
  for (i in 1:length(trans)){ # Find degeneracy in every position
    mpos <- ceiling(i/3)
    ncod <- as.numeric(PopGenome:::codonise64(a[mpos,,drop=F]))
    # test <- c(test,code[ncod])
    fold <- c(fold,assocmat[ncod,indices[i]])
    #if (tab$strand[1] == "-") {fold <- rev(fold)}
  }
  pos1 <- 1

  # Checking fold  
  
  for (n in 1:nrow(tab)){
    pos2 <- pos1+tab[n,2]-tab[n,1] # Position in the fold object (protein)
    local.fold <- fold[pos1:pos2]
    gffseq[tab[n,1]:tab[n,2]][local.fold == 0] <- 0
    # Four-fold are less constrained that 5UTR, 3UTR and introns:
    constrained <- gffseq[tab[n,1]:tab[n,2]] == 5 | gffseq[tab[n,1]:tab[n,2]] == 3 | gffseq[tab[n,1]:tab[n,2]] == 9
    gffseq[tab[n,1]:tab[n,2]][!UTR & local.fold == 4] <- 4
    pos1 <- pos2
  }
  print(x)
}
Sys.time()-init

#### ALTERNATIVE METHOD:

# Account for the forward and reverse strands:
START <- cdsgff[,3]
END <- cdsgff[,4]
REV <- cdsgff[,5]=="-"
pos.strand.shift <- as.integer(cdsgff$frame[!REV]) # 3rd column is the shift (e.g. 2) from reading frame
neg.strand.shift <- as.integer(cdsgff$frame[REV])

START[!REV] <- START[!REV] + pos.strand.shift # We add the shift to the START coordinates in non-reverse regions
END[!REV] <- END[!REV] + as.numeric(gsub(3,0,3-pos.strand.shift)) # Subtract from reverse when ORF is not 0

START[REV] <- cdsgff[,4][REV] - rev.strand.shift # We take the END coordinates for reverse regions and subtract shift
END[REV] <- cdsgff[,3][REV] - as.numeric(gsub(3,0,3-neg.strand.shift)) 

# Extract the positions corresponding to the CDS
CDS <- unique(data.frame(START,END,REV))
all.codons <-apply(CDS,1,function(x){CHR[x[1]:x[2]]}) # Select positions specified in codons from FASTA file 
altenv <- new.env()
altenv$count <- 1 # Counter for lapply
all.codons <- lapply(all.codons,function(x){
  if (CDS[altenv$count,3] == 1){
    x <- komplement[x]
  }
  altenv$count <- altenv$count+1
  return(x)
}) # We obtain a list containing the numeric code for every codon

s <- lapply(all.codons,function(x){
  fold <- vector() # Contains the degeneracy fold in each position
  test <- vector() # Contains translated aminoacids
  x2 <- x[1:(length(x)-length(x)%%3)] # Sequences that are not multiple of 3
  indices <- seq(0,length(x2)-1)%%3+1 # +1 because we are using assocmat
  a <- matrix(x2,ncol=3,byrow=TRUE)
  for (i in 1:length(x2)){
    mpos <- ceiling(i/3)
    ncod <- as.numeric(PopGenome:::codonise64(a[mpos,,drop=F]))
    test <- c(test,code[ncod])
    fold <- c(fold,assocmat[ncod,indices[i]])
  }
  return(fold)
})

for (i in 1:length(feature)) {
  gffseq[start(feature[i]):end(feature[i])] <- symbol
}

###########################
## GENE-BY-GENE ANALYSIS ##
###########################

gff.table[gff.table$feature == "gene",]



##############
## CHECKING ##
##############

gene.names <- gff.table[gff.table$feature == "gene",]$ID
cds.names <- sapply(unique(cdsgff$ID),function(x){substr(x,1,nchar(x)-2)})
sum(unique(cds.names)) > sum(unique(gene.names))
gene.names %in% cds.names # Not all genes have CDS... Some are not translated? 
# The most plausible explanation is that they are exons from non-coding RNAs

# Start and end of txDB and gff are the same. Therefore, the differences between 
# the number of nucleotides are due to something else
start(genes) %in% gff.table[gff.table[,2] == "gene",3]
end(genes) %in% gff.table[gff.table[,2] == "gene",4]
