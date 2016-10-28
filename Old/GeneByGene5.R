# setwd("~/Documents/2_GenomicsData/TestGeneByGene")
## filename <- "merge.vcf.gz"; ini <- 38944867; end <- 38966701; chrom <- "22"
# gpfile <- "chr22_filtered.vcf.gz" ; alnfile <- "chr22_aln.vcf.gz"

# UPDATE: Instead of segmenting the PopGenome object, we will remove those positions that
# are missing via vector.
# UPDATE 2: 

library(PopGenome)
library(GenomicRanges)
library(stringr)

##########################
## FUNCTION DECLARATION ##
##########################

## POLYMORPHISM FUNCTIONS ##

Pi <- function(k,m,n) { # Window number
  comb <- choose(n,2) # Binomial coefficient = combination without repetition
  pi <- k/(comb*m)
  return(round(pi,7))
}
Theta <- function(S,m,n) { # S=Segregating sites; m=total sites; n=number of samples
  summ <- 0
  for (i in 1:(n*2-1)) {
    summ <- summ + 1/i
  }
  theta <- (S/m)/summ
  return(round(theta,7))
}

## MKT CALCULATION ##

mkt.extended <- function (sel=0,neu=4) {
  
  # STANDARD: Same result as MKT function if done with all syn.
  Pneu <- sum(polym[bial.class==neu],na.rm=TRUE) # For all syn, >=1
  Psel <- sum(polym[bial.class==sel],na.rm=TRUE)
  Dneu <- sum(divtotal[bial.class==neu],na.rm=TRUE)
  Dsel <- sum(divtotal[bial.class==sel],na.rm=TRUE)
  NI <- (Psel/Pneu)/(Dsel/Dneu)
  alpha <- 1-NI
  DoS <- Dsel/(Dsel+Dneu)-Psel/(Psel+Pneu) # Direction of selection
  contingency.std <- matrix(c(Psel,Pneu,Dsel,Dneu),c(2,2))
  test <- if(!is.na(sum(contingency.std))){
    if(sum(contingency.std)>0){fisher.test(contingency.std)$p.value}
  } else {NA}
  
  # ADAPTED:
  
  # 1) Use four-fold degenerated sites 2) Remove sites according to DAF
  
  # Using DAF results in a loss of information, as it is not available for some
  # variants. For MAF, simply retrieve from MAF.
  
  neuMAF <- as.vector(na.omit(MAF[bial.class==neu & polym==TRUE]))
  selMAF <- as.vector(na.omit(MAF[bial.class==sel & polym==TRUE]))
  
  Pneu.less5 <- sum(neuMAF<0.05) # P0 MAF less than 5%
  Pneu.more5 <- sum(neuMAF>0.05) # P0 MAF more than 5% 
  Psel.less5 <- sum(selMAF<0.05) # Pi MAF less than 5% 
  Psel.more5 <- sum(selMAF>0.05) # Pi MAF more than 5%
  
  Psel.neutral.less5 <- Psel*(Pneu.less5/Pneu) # Proportion of neutral within the MAF < 5% class
  Psel.neutral <- Psel.neutral.less5 + Psel.more5 # For alpha
  Psel.weak <- Psel.less5 - Psel.neutral.less5
  
  alpha.cor <- 1-(Psel.neutral/Pneu)*(Dneu/Dsel)
  contingency <- matrix(c(Psel.neutral,Pneu,Dsel,Dneu),c(2,2))
  test.cor <- if(!is.na(sum(contingency))){
    if(sum(contingency)>0){fisher.test(contingency)$p.value}
  } else {NA}
  
  # OTHER ESTIMATORS:
  
  # To fully implement extended MKT, we need ms/mns, that is, the number of sites
  # of each class. Doing that would require modifying the 'set.synnonsyn2' function
  # to obtain all codons and check fold in every position (0,1,2)
  
  m.neu <- sum(bial.class == neu,na.rm=T)
  m.sel <- sum(bial.class == sel,na.rm=T)
  
  f <- (m.neu*Psel.neutral)/(m.sel*Pneu) # Neutral sites
  b <- (Psel.weak/Pneu)*(m.neu/m.sel)
  y <- (Psel/Pneu-Dsel/Dneu)*(m.neu/m.sel)
  d <- 1 - (f+b)
  
  return(c(Psel,Pneu,Dsel,Dneu,alpha,test,Psel.neutral,alpha.cor,test.cor,DoS,f,b,y,d))
}

########################
## POPGENOME ANALYSIS ## 
########################

popanalysis <- function(filename,ini,end)
  region <- readVCF(filename,numcols=5000,tid=chrom,from=ini,to=end,include.unknown=TRUE,gffpath=sprintf("chr%s.gff",chrom))
  # region <- set.synnonsyn(region,ref.chr=sprintf("chr%s.fa",chrom),save.codons=FALSE) 
  # Syn-nonsyn is not needed if we only use 0- and 4-fold
  
  ## DEFINE POPULATION/OUTGROUP ##
  
  # The outgroup can be defined as a different population so that calculations of diversity and neutrality 
  # are performed separately. 
  
  individuals <- get.individuals(region)[[1]]
  humans <- individuals[which(individuals != c("Chimp","Chimp.2"))]
  chimp <- individuals[which(individuals == c("Chimp","Chimp.2"))]
  region <- set.populations(region,list(humans,chimp))
  region <- set.outgroup(region,chimp)
  
  ## FILTER THE BIALLELIC MATRIX ##
  
  n <- length(humans) # N = samples. Diploid, not divided by 2
  total <- length(individuals) # Total number of samples (including outgroup)
  bial <- get.biallelic.matrix(region,1) # Biallelic matrix of the window
  if (is.null(bial)||dim(bial)[2]==0) { # When no variants are detected
    newrow <- rep(0,13) # Empty rows
    return(newrow)
  }
  ac.bial <- as.numeric(colnames(bial)) %in% ac.pos # Accessible positions in biallelic matrix
  bial <- bial[,ac.bial]
  wsize <- end-ini+1 # GFF coordinates, 1-based
  
  # Bial contains only positions that are found in the accessibility mask
  
  bialhuman <- bial[1:n,,drop=F] # Remove outgroup (drop = F to keep 1 dimension)
  misshuman <- colSums(is.na(bialhuman))>0 # Sites missing in humans (e.g. structural variants)
  polym <- apply(bialhuman,2,sum)>0 & !misshuman # Sites polymorphic in humans w/o missing
  bialhuman <- bialhuman[,polym,drop=FALSE] # Keep only polymorphic sites
  
  ## SITE FREQUENCY SPECTRUM (FOLDED) ##
  
  MAF <- colSums(bial[1:n,,drop=F] == 1)/nrow(bialhuman) # We DO NOT remove monomorphic alleles
  MAF.df <- data.frame(POS=names(MAF),MAF=unname(MAF)) # Position and DAF of biallelic variants
  
  # IDENTIFY ANCESTRAL ALLELES (AA) FOR SNPS IN THE REGION
  temp <- system(sprintf('gunzip -c %s | cut -f2,4,8',filename),intern=TRUE) # Extract all variants from the merge.vcf.gz file  
  gpimport <- read.table(textConnection(temp),sep="\t",stringsAsFactors = FALSE,header=TRUE)
  # gpimport <- read.table("gpimport.txt",sep="\t",stringsAsFactors = FALSE,header=TRUE)
  gpances <- toupper(sapply(gpimport[,3],function(x){str_match(x,"AA=([:alpha:])\\|")[2]},USE.NAMES=FALSE)) # Capture the ancestral allele
  aarefs <- na.omit(cbind(gpimport[,1:2],ANC=gpances)) # Na.omit removes positions without known AA and divergent sites
  
  # CALCULATE DERIVED ALLELE FREQUENCY (DAF)
  mafan <- merge(MAF.df,aarefs) # Intersection between SNPs in the biallelic matrix and the AA in the VCF file
  if (dim(mafan)[1] > 0) { # Avoid windows where variants have no known AA
    mafan$DAF <- NA # DAF is equal to MAF if derived = minor, but 1 - MAF if otherwise
    mafan[mafan$REF != mafan$ANC,]$DAF <- 1- mafan[mafan$REF != mafan$ANC,]$MAF 
    mafan[mafan$REF == mafan$ANC,]$DAF <- mafan[mafan$REF == mafan$ANC,]$MAF 
  }
  mafan <- merge(MAF.df,mafan,all.x=TRUE) # Retain originals
  
  winDAF <- as.vector(na.omit(mafan$DAF[windex]))
  DAF <- hist(winDAF,seq(0.0,1,0.05),plot=F)$counts
  DAF <- paste(DAF,collapse = ";")
  
  ## POLYMORPHISM ## 
  
  # DETERMINE M (EXCLUDING MISSING AND POLYALLELIC)
  # Total number of sites: 1) Remove NA in humans 2) Remove polyallelic sites
  if (!length(region@region.data@polyallelic.sites) == 0) { # Make sure list exists
    polysites <- sum(region@region.data@polyallelic.sites[[1]]) # Positions of all polyalleles
  } else {polysites <- 0} # If not available, we assume 0
  m <- wsize-sum(misshuman,na.rm=T)-polysites
  # DETERMINE S AND K (WHEN VARIANTS ARE AVAILABLE)
  if (is.null(bialhuman)||dim(bialhuman)[2]==0) {
    S <- 0
    k <- 0
  } else {
    S <- ncol(bialhuman) # Number of variants
    freqs <- apply(bialhuman,2,table)
    k <- sum(freqs[1,]*freqs[2,]) # Note that k and K are different!
  }
  
  # Functions are declared above!
  
  ## DIVERGENCE ##
  
  # COUNT THE NUMBER OF UNKNOWNS IN THE OUTGROUP:
  unknowns <- sum(is.nan(bial[total,])) # Displayed as NaN in the biallelic matrix
  
  # COUNT THE NUMBER OF DIVERGENCE SITES:
  divtotal <- bial[total,] == 1 & !polym # 1/1 in the outgroup and excluding polymorphic (cannot be counted for divergence)
  divsites <- sum (divtotal,na.rm=T) # Number of divergent sites
  # CALCULATE D AND K:
  mout <- m - unknowns
  D <- round(divsites/mout,7) # Observed divergence: Proportion of sites with divergent nucleotides
  K <- round(-3/4*log(1-4/3*D),7) # Real divergence: Jukes and Cantor model to account for multiple hits
  
  ## NATURAL SELECTION REGIMES ##
  
  sites <- as.numeric(colnames(bial)) # Biallelic sites
  bial.class <- gffseq[sites] # GFF feature in each site. Equivalent to the former 'syn'
  # syn <- region@region.data@synonymous[[1]][ac.bial] # Contains syn and non-syn positions [accessible]
  
  # STORING MKT RESULTS FOR DIFFERENT FUNCTIONAL CLASSES:
  
  mkt.0fold.4fold <- mkt.extended(sel=0,neu=4) # sel = 0-fold; neu= 4-fold
  mkt.introns.4fold <-  mkt.extended(sel=9,neu=4) # sel = exon; neu= 4-fold
  mkt.5UTR.4fold <-  mkt.extended(sel=5,neu=4) # sel = 5-UTR; neu= 4-fold
  mkt.3UTR.4fold <-  mkt.extended(sel=3,neu=4) # sel = 3-UTR; neu= 4-fold
  bial.class[bial.class == 3] <- 5 # To combine all UTR
  mkt.UTR.4fold <-  mkt.extended(sel=5,neu=4) # sel = exon; neu= 4-fold
  bial.class[is.na(bial.class)] <- 1 # New code for all intergenic (instead of NA)
  mkt.inter.4fold <-  mkt.extended(sel=1,neu=4) # sel = exon; neu= 4-fold
  
  ## ADD NEW ROW ##
  newrow <- c(S,Pi(k,m,n),DAF,divsites,D,K,unknowns,mkt.0fold.4fold,mkt.introns.4fold,mkt.5UTR.4fold,mkt.3UTR.4fold,mkt.UTR.4fold,mkt.inter.4fold )
  
  return(tabsum)  
}

###################
## GENE ANALYSIS ##
###################

## RETRIEVE ENTREZ GENES FROM TXDB:

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)

load("gffseq_chr22.RData")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb) <- sprintf("chr%s",chrom)
grgenes <- genes(txdb)
start(grgenes) <- start(grgenes)-500 # Upstream(-strand)
end(grgenes) <- end(grgenes)+500 # Downstream(+strand)

# Add +500, -500 

## CREATE A GENE SPECIFIC MERGE:

merge.vcf <- function(ini,end) {
  system(sprintf("bcftools merge -Oz --missing-to-ref -o merge_gene.vcf.gz -r %s:%d-%d %s %s",
                 chrom,ini,end,gpfile,alnfile))
  system("tabix -p vcf merge_gene.vcf.gz")
  filename <- "merge_gene.vcf.gz"
  return(filename) 
}

# Note that the coordinates of the MASK are in BED format and therefore 0-based, whereas the GFF with
# the genes and the GRanges objects are 1-based. To convert from 0 to 1-based, START+1:END.

# TABLE CONTAINING THE DATA:
tabsum <- as.data.frame(matrix(numeric(length(grgenes)*13),ncol=13,nrow=length(grgenes)))
colnames(tabsum) <- c("S","Pi","DAF","Divsites","D","K","Unknown","Alpha","Fisher","Pns","Ps","Dns","Ds")

## USING THE FASTA SEQUENCE:
library("Biostrings")
maskfasta <- readBStringSet("chr22.pilot_mask.fasta.gz") # Reading files in gz format IS supported
init <- Sys.time()
missing <- numeric(length(grgenes))
for(i in 1:length(grgenes)) {
  ini <- start(grgenes[i]); end <- end(grgenes[i])
  mask.local <- strsplit(as.character(subseq(maskfasta,start=ini,end=end)),"")[[1]]
  pass <- mask.local == "P"
  if (sum(pass) == 0) {missing[i] <- 100}
  ac.pos <- (ini:end)[pass] # Vector with gene positions that are accessible
  missing[i] <- (1-sum(pass)/length(pass))*100 # Proportion of positions that do not                                                                                                                                                                                                                   pass the filter
  filename <- merge.vcf(ini,end)                                                                                                                        
}
Sys.time()-init


###################
## FINAL TOUCHES ##
###################

library(biomaRt)
mart <- useMart(host="www.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL")
atts <- c("entrezgene","hgnc_symbol") # Retrieve the HUGO-approved gene symbol
gene.names <- genes$gene_id
annotation.h <- getBM(attributes=atts,uniqueRows=FALSE,filters="entrezgene",values=gene.names,mart=mart)
