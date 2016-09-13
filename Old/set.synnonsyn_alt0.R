# ref.chr = "chr22.fa"; xyz = 1; object = region; save.codons = FALSE; four.fold = TRUE

setGeneric("set.synnonsyn2", function(object,ref.chr=FALSE,save.codons=FALSE,four.fold=TRUE) standardGeneric("set.synnonsyn2"))
 setMethod("set.synnonsyn2", "GENOME",
 function(object,ref.chr,save.codons,four.fold){

codonise64 <- function(s){ # Function to identify the codon number (by BPfeifer). Should be integrated in PopGenome package.
 
 m<-dim(s)[2]
 n<-dim(s)[1]
 
 if((m%%3)>0){
   stop("length of coding sequence cannot divide by 3!") 
 }
 
 t <- matrix(,n,m/3)
 p <- seq(1,m,3)
 q <- 1:(m/3)
 
 s[s>=5] <- 100 # if gaps or unknown in the codon
 
 t[,q] <- (s[ ,p]-1)*16 + (s[ ,p+1]-1)*4 + (s[ ,p+2]-1)+1
 t[t>64]<-65
 
 return(t)
}   
   
# codontable() # 1st matrix - aa of the triplets, each row is a different code (first standard)
# 2nd matrix - nucleotide combinations of each aa (row), matching the matrix above. U=1, C=2, A=3, G=4.
# region@region.data@CodingSNPS - Coding SNPs TRUE/FALSE

check_fold <- function(triplet,pos) {
 # CHECK THAT NO ARGUMENTS ARE MISSING
 if (missing(triplet) | missing(pos)){
   stop("Arguments are missing")
 }
 # EVALUATE INPUT VARIABLES:
 T <- codontable()
 trips <- T$Triplets
 code <- T$Protein[1,]
 comb <- c(codonise64(triplet)) # Assumes matrix 1x3. Vectorize it with c()
 aa <- code[comb] # Original aminoacid of the triplet
 # EVALUATE DEGENERACY:
 fold <- 0
 for (i in 1:4) {
   triplet[pos+1] <- i # Pos in cod.pos is 0/1/2, but indices in R start at 1
   comb <- c(codonise64(triplet)) # Assumes matrix 1x3. Vectorize it with c()
   altaa <- code[comb]
   # cat(c(triplet,",",altaa,"\n"))
   if (altaa == aa) {fold <- fold+1}
 }
  return (fold == 4)
}
   
if(ref.chr[1]==FALSE){
stop("Please verify the reference sequence")
}

if(object@gff.info == FALSE) {
stop("No GFF file was read in !")
}


for (xyz in 1:length(ref.chr)){
# erstmal nur fuer ein Chunk

Coding.matrix            <- object@region.data@Coding.matrix2[[xyz]][,] # weil ff object, 2 because (fitting GFF)
biallelic.sites2         <- object@region.data@biallelic.sites2[[xyz]]  # Number of variant in the VCF file. E.g. 2 - variant n 2.
biallelic.sites          <- object@region.data@biallelic.sites[[xyz]] # Genomic coordinates of each variant
START                    <- object@region.data@reading.frame[[xyz]][,] #  Coords of coding region and reading frame
REV                      <- object@region.data@rev.strand[[xyz]][,]    #  reverse strand information for coding regions
REV                      <- as.logical(REV)  # If REV = TRUE, that coding region is in the reverse strand
CodingSNPS               <- object@region.data@CodingSNPS[[xyz]]

if(sum(CodingSNPS)==0){
warning("No coding SNPs in this region !")
print(object@region.names[xyz])
next
}

# The START position is altered depending whether coding regions are reverse

#print(REV)
if(any(!REV)){
pos.strand.shift  <- START[!REV,3] # 3rd column is the shift (e.g. 2) from reading frame

#zero			 <- pos.strand.shift==0
#one 			 <- pos.strand.shift==1
#two 			 <- pos.strand.shift==2	
#pos.strand.shift[zero]   <- 0
#pos.strand.shift[one]    <- 2
#pos.strand.shift[two]    <- 1

}else{
pos.strand.shift         <- 0
}

if(any(REV)){ 
rev.strand.shift         <- START[REV,3] # 3rd column is the shift (e.g. 2) from reading frame

#zero			 <- rev.strand.shift==0 
#one 			 <- rev.strand.shift==1 
#two 			 <- rev.strand.shift==2
#rev.strand.shift[zero]   <-  0
#rev.strand.shift[one]    <-  2
#rev.strand.shift[two]    <-  1
}else{
rev.strand.shift         <- 0
}

START                    <- START[,1] # Starting positions of the coding regions
START[!REV]              <- START[!REV] + pos.strand.shift # We add the shift to the START coordinates in non-reverse regions

START[REV]               <- object@region.data@reading.frame[[xyz]][REV,2] - rev.strand.shift # We take the END coordinates for reverse regions and subtract shift

Coding.matrix            <- Coding.matrix # Relative start/end of coding regions

# define an evironment
synGLOBAL <- new.env() 

# Create Region and save size of region 
synGLOBAL$SIZE  <- numeric(dim(Coding.matrix)[1])
synGLOBAL$count <- 1

erg  <- apply(Coding.matrix,1,function(xx){ # List containing indices of coding regions

 region                            <- xx[1]:xx[2] # Coding positions per region
 synGLOBAL$SIZE[synGLOBAL$count]   <- length(region) # Length of each region
 synGLOBAL$count                   <- synGLOBAL$count + 1 # Number of regions
 return(region)

})

# What are the real positions? Genomic coordinates of coding SNPs.
 erg         <- unlist(erg)

 bial.pos    <- match(erg,biallelic.sites2) #.Call("my_match_C",erg,biallelic.sites2)
 bial.pos[bial.pos==-1] <- NaN
 #return(bial.pos)
 bial.pos    <- biallelic.sites[bial.pos]
 
# RogerNote: Isn't it possible to just do? 
# bial.pos <- biallelic.sites[CodingSNPS]
# No! Missing NA positions of non-biallelic 

#print(bial.pos)
#print(biallelic.sites2)

# Create Start Vector
 synGLOBAL$count <- 1
 vec <- sapply(START,function(x){
      gg              <- rep(x,synGLOBAL$SIZE[synGLOBAL$count])       
      synGLOBAL$count <- synGLOBAL$count + 1

 return(gg)
 })

# Create REV Vector #----
 synGLOBAL$count <- 1
 vec_rev <- sapply(REV,function(x){
      gg              <- rep(x,synGLOBAL$SIZE[synGLOBAL$count])       
      synGLOBAL$count <- synGLOBAL$count + 1

 return(gg)
 })


 START.vec   <- unlist(vec)
 REV.vec     <- unlist(vec_rev) #---

#print(length(REV.vec))
#print(length(START.vec))
#print(length(bial.pos))

 # DEFINE CODONS #
 
 # START contains the starting position of the coding region, so any other
 # position can be referenced to it. With modulo we find the remainder of the
 # quotient, which indicates the position in a triplet: 0 is first, 1 is second...
 
 cod.pos <- (bial.pos - START.vec)%%3 
 # in case of reverse strand 
 cod.pos[REV.vec] <- (START.vec[REV.vec]-bial.pos[REV.vec])%%3


 # Delete NaNs 
 cod.pos     <- cod.pos[!is.na(bial.pos)]
 rev.pos     <- REV.vec[!is.na(bial.pos)] #----
 bial.pos    <- bial.pos[!is.na(bial.pos)]

 # Now we have the position of each variant in the triplet!
 
 #print(START.vec[300])

# print(length(bial.pos))
# print(length(cod.pos))
# print(length(rev.pos))

 ids         <- !duplicated(bial.pos)
 cod.pos     <- cod.pos[ids]
 bial.pos    <- bial.pos[ids]
 rev.pos     <- rev.pos[ids] #----

# print(bial.pos)
# print(length(cod.pos))
# print(length(rev.pos))

 # print(rev.pos[1215])
 # print(bial.pos[1215])
 # print(cod.pos[1215])
 #print(rev.pos)

# CREATE THE CODONS
# bial.pos and cod.pos: create triplets around bial.pos

codons <- matrix(,length(cod.pos),3) # As many rows as variants; 3 columns for triplets

for (xx in 1:length(cod.pos)){ # Assigns coordinates to each codon
   
    if(rev.pos[xx]){# reverse strand #FIXME
     if(cod.pos[xx]==0){codons[xx,]=c(bial.pos[xx],bial.pos[xx]-1,bial.pos[xx]-2);next}
     if(cod.pos[xx]==1){codons[xx,]=c(bial.pos[xx]+1,bial.pos[xx],bial.pos[xx]-1);next}
     if(cod.pos[xx]==2){codons[xx,]=c(bial.pos[xx]+2,bial.pos[xx]+1,bial.pos[xx]);next}
    }else{
     if(cod.pos[xx]==0){codons[xx,]=c(bial.pos[xx],bial.pos[xx]+1,bial.pos[xx]+2);next}
     if(cod.pos[xx]==1){codons[xx,]=c(bial.pos[xx]-1,bial.pos[xx],bial.pos[xx]+1);next}
     if(cod.pos[xx]==2){codons[xx,]=c(bial.pos[xx]-2,bial.pos[xx]-1,bial.pos[xx]);next}
     }

}
#print(cod.pos)
#print(codons)

## Reading the reference chromosome
file.info <- .Call("get_dim_fasta",ref.chr[xyz]) # .Call invokes the C function get_ind_fasta (src folder): start and end positions of chr.

gc()
#print(file.info)

CHR       <- .Call("get_ind_fasta",ref.chr,1,file.info[[1]][2]) # .Call invokes get_ind_fasta, which extracts the letters in FASTA to numbers.

#print(CHR[1:10])

# Create codons with nucleotides
Nuc.codons    <- CHR[codons] # Select positions specified in codons from FASTA file
Nuc.codons    <- matrix(Nuc.codons,ncol=3) # Put them in a 3-col matrix

#print(Nuc.codons)

ALT           <- Nuc.codons
REF           <- Nuc.codons
Subst         <- object@region.data@biallelic.substitutions[[xyz]] # All substitions in the matrix
minor         <- Subst[1,CodingSNPS]
major         <- Subst[2,CodingSNPS]

komplement <- c(4,3,2,1,5) # For regions in the reverse strand: opposite of T=1,C=2,G=3,A=4

ffold <- rep(NA,dim(Nuc.codons)[1]) # To check whether it is 4fold degenerate

for(xx in 1: dim(Nuc.codons)[1]){ # For each codon with variants
 if(rev.pos[xx]){

  # Convert to komplement nucleotides
  REF[xx,]  <- komplement[REF[xx,]]
  ALT[xx,]  <- komplement[ALT[xx,]]
  minor[xx] <- komplement[minor[xx]]
  major[xx] <- komplement[major[xx]]
  ###########
  
  # Why does he use the minor allele for REF? Then the codon change is backwards.
  # Ask Sonia. It is done on purpose, the getcodons function indicates it correctly.
  if(cod.pos[xx]==0){REF[xx,1] <- minor[xx];ALT[xx,1]<-major[xx]}
  else if(cod.pos[xx]==1){REF[xx,2] <- minor[xx];ALT[xx,2]<-major[xx]}
  else if(cod.pos[xx]==2){REF[xx,3] <- minor[xx];ALT[xx,3]<-major[xx]}
 }else{
  if(cod.pos[xx]==0){REF[xx,1] <- minor[xx];ALT[xx,1]<-major[xx]}
  else if(cod.pos[xx]==1){REF[xx,2] <- minor[xx];ALT[xx,2]<-major[xx]}
  else if(cod.pos[xx]==2){REF[xx,3] <- minor[xx];ALT[xx,3]<-major[xx]}
 }
  
  # CHECK 4-FOLD DEGENERATE #

  pos <- cod.pos[xx]
  ffold[xx] <- check_fold(REF[xx,,drop=F],pos) # We use a matrix of 1 dimension

} 

#print(REF)
#print(ALT)


# Coding Codons ...

ALT <- codonise64(ALT)
REF <- codonise64(REF)

if(save.codons){
saveALTREF <- cbind(REF,ALT)
}

CC  <- codontable()

ALT <- CC$Protein[1,ALT]
REF <- CC$Protein[1,REF]

# Check differences between REF and ALT codons

CHECK <- cbind(ALT,REF)

erg <- apply(CHECK,1,function(x){return(length(unique(x)))}) # Length of unique vector: 2 if aa changes
erg[erg==2] <- 0 #nonsyn
# Label the 4fold degenerate sites instead of synonymous:
if (four.fold){
  erg[erg==1 & ffold] <- 4 # Four-fold sites as 4
  erg[erg==1 & !ffold] <- 1 # Synonym sites as 1
  
} else {
  erg[erg==1] <- 1 #syn
  }

# SAVING CHANGES TO OBJECT GENOME #

# Change object of class GENOME
change <- object@region.data
change@synonymous[[xyz]][CodingSNPS] <- erg

### save codons
if(save.codons){
n.coding.snps     <- sum(CodingSNPS)
codonlist  <- vector("list",n.coding.snps)
count <- 1

for(vv in 1:n.coding.snps){
    codonlist[[vv]] <- saveALTREF[vv,]
}
change@codons[[xyz]] <- codonlist
}
#######################

object@region.data <- change 
}# End Iteration over chunks or chromosomes

return(object)

})
