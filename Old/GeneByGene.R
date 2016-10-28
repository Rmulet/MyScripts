# setwd("~/Documents/2_GenomicsData/TestGeneByGene")
# filename <- "merge.vcf.gz"; ini <- 38944867; end <- 38966701; chrom <- "22"
# gpfile <- "chr22_filtered.vcf.gz" ; alnfile <- "chr22_aln.vcf.gz"

# Put the 1000GP-ALN merge in a temp???
library(PopGenome)

load("gffseq_chr22.RData")

## CREATE A GENE SPECIFIC MERGE:

system(sprintf("bcftools merge -Oz --missing-to-ref -o merge_gene.vcf.gz -r %s:%d-%d %s %s",
               chrom,ini,end,gpfile,alnfile))
system("tabix -p vcf merge_gene.vcf.gz")
filename <- "merge_gene.vcf.gz"

# Note that the coordinates of the MASK are in BED format and therefore 0-based, whereas the GFF with
# the genes and the GRanges objects are 1-based. To convert from 0 to 1-based, START+1:END.

## USING GRANGES:
library(GenomicRanges)

genes <- gff.table[gff.table$feature == "gene",c(1:5,7)]
genes[,3:4] <- c(genes$start - 500,genes$end + 500) # Expand to up- and down-stream regions
grgenes <- with(genes,GRanges(chr,IRanges(start,end)))
mask <- read.table("pilot_mask.chr22.bed",colClasses = c("character","integer","integer","NULL"))
colnames(mask) <- c("chr","start","end")
mask$chr <- sapply(strsplit(mask$chr,"chr"),"[[",2)
mask$start <- mask$start+1 # Add +1 to convert from BED to GFF
grmask <- with(mask,GRanges(chr,IRanges(start,end)))

overlap <- countOverlaps(query=grgenes,subject=grmask)
grgenes <- grgenes[overlap > 0] # Genes that are at least partially found in the accessible regions

init <- Sys.time()
for(i in 1:length(grgenes)) {
  inter <- intersect(grgenes[i],grmask)
  if (sum(width(inter)) < 0.5*width(grgenes[i])) {
    print (sum(width(inter)))
  }
  ini <- start(grgenes[1]); end <- end(grgenes[1])}
  
Sys.time()-init

seq <- gffseq[ini:end]

region <- readVCF(filename,numcols=5000,tid=chrom,from=ini,to=end,include.unknown=TRUE,gffpath=sprintf("chr%s.gff",chrom))
region <- set.synnonsyn(region,ref.chr=sprintf("chr%s.fa",chrom),save.codons=FALSE) 

## USING THE FASTA SEQUENCE:
library("Biostrings")
maskfasta <- readBStringSet("chr22.pilot_mask.fasta")
for(i in 1:nrow(genes)) {
  mask.local <- strsplit(as.character(subseq(f,start=genes[1,3],end=genes[1,4])),"")[[1]]
  pass <- mask.local == "P"
  if (sum(pass) < 50) # Threshold of positions that pass the filter
  
  if (sum(width(inter)) < 0.5*width(grgenes[i])) {
    print (sum(width(inter)))
  }
  ini <- start(grgenes[1]); end <- end(grgenes[1])}


qgene <- m[start(grgenes[1]):end(grgenes[1])] # Gene quality




########################
## POPGENOME ANALYSIS ## 
########################

# We can either concatenate regions or calculate polymorphism manually

## DEFINE POPULATION/OUTGROUP ##

# The outgroup can be defined as a different population so that calculations of diversity and neutrality 
# are performed separately. 

load("gffseq_chr22.RData")

individuals <- get.individuals(region)[[1]]
humans <- individuals[which(individuals != c("Chimp","Chimp.2"))]
chimp <- individuals[which(individuals == c("Chimp","Chimp.2"))]
region <- set.populations(region,list(humans,chimp))
region <- set.outgroup(region,chimp)

## NEUTRALITY TESTS ##

# This calculates neutrality stats (Watterson's Theta, Fu and Li, Tajima's D) and the
# number of segregating sites (S). However, S excludes missing information in the outgroup.
region <- neutrality.stats(region)

## IDENTIFY BIALLELIC MATRIX ## 

n <- length(humans) # N = samples. Diploid, not divided by 2
total <- length(individuals) # Total number of samples (including outgroup)
bial <- get.biallelic.matrix(region,1) # Biallelic matrix of the window
if (is.null(bial)||dim(bial)[2]==0) { # When no variants are detected
  newrow <- rep(0,13) # Empty rows
  tabsum[window,] <- newrow
  next
}
bialhuman <- bial[1:n,,drop=F] # Remove outgroup (drop = F to keep 1 dimension)
misshuman <- colSums(is.na(bialhuman))>0 # Sites missing in humans (e.g. structural variants)
polym <- apply(bialhuman,2,sum)>0 & !misshuman # Sites polymorphic in humans w/o missing
bialhuman <- bialhuman[,polym,drop=FALSE] # Keep only polymorphic sites

# TABLE CONTAINING THE DATA:
tabsum <- as.data.frame(matrix(numeric(nwin*13),ncol=13,nrow=nwin))
colnames(tabsum) <- c("S","Pi","DAF","Divsites","D","K","Unknown","Alpha","Fisher","Pns","Ps","Dns","Ds")

# DETERMINE SEGREGATING SITES (EXCLUDING OUTGROUP)

# DETERMINE M (EXCLUDING MISSING AND POLYALLELIC)
winstart <- windows[window,1] # Select the start position in that window
# Total number of sites: 1) Remove NA in humans 2) Remove polyallelic sites
if (!length(region@region.data@polyallelic.sites) == 0) { # Make sure list exists
  polyal <- region@region.data@polyallelic.sites[[1]] # Positions of all polyalleles
  polysites <- sum(!is.na(match(polyal,winstart:(winstart+wsize)))) # N polyallelic sites in that window
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

## POLYMORPHISM ## 
# PopGenome has its own 'neutrality' and 'diversity' stats modules, but diversity takes a 
# very long time (55 seconds including Pi for a region of 100 kb). 
Pi <- function(k,m,n) { # Window number
  comb <- choose(n,2) # Binomial coefficient = combination without repetition
  pi <- k/(comb*m)
  return(round(pi,7))
}

## DIVERGENCE ##

# COUNT THE NUMBER OF UNKNOWNS IN THE OUTGROUP:
unknowns <- sum(is.nan(bial[total,])) # Displayed as NaN in the biallelic matrix

# COUNT THE NUMBER OF DIVERGENCE SITES:
divtotal <- bial[total,] == 1 & !polym # 1/1 in the outgroup and excluding polymorphic
divsites <- sum (divtotal,na.rm=T) # Polymorphic cannot be counted for divergence
# CALCULATE D AND K:
mout <- m - unknowns
D <- round(divsites/mout,7) # Observed divergence: Proportion of sites with divergent nucleotides
K <- round(-3/4*log(1-4/3*D),7) # Real divergence: Jukes and Cantor model

## NATURAL SELECTION REGIMES ##

winsites <- allwinsites[[window]] # Biallelic sites in this window
bial.gff <- gffseq[allwinsites[[window]]] # GFF feature in each site
syn <- slide@region.data@synonymous[[window]] # Contains syn and non-syn positions
syn[bial.gff == 4] <- 4

## MKT CALCULATION

# STANDARD: Same result as MKT function if done with all syn.

Ps <- sum(polym[syn==4],na.rm=TRUE) # For all syn, >=1
Pns <- sum(polym[syn==0],na.rm=TRUE)
Ds <- sum(divtotal[syn==4],na.rm=TRUE)
Dns <- sum(divtotal[syn==0],na.rm=TRUE)
NI <- (Pns/Ps)/(Dns/Ds)
alpha <- 1-NI
DoS <- Dns/(Dns+Ds)-Pns/(Pns+Ps) # Direction of selection
contingency.std <- matrix(c(Pns,Ps,Dns,Ds),c(2,2))
test <- if(!is.na(sum(contingency.std))){
  if(sum(contingency.std)>0){fisher.test(contingency.std)$p.value}
} else {NA}

# ADAPTED:

# 1) Use four-fold degenerated sites 2) Remove sites according to DAF

# Using DAF results in a loss of information, as it is not available for some
# variants. For MAF, simply retrieve from freqs.df.

synDAF <- as.vector(na.omit(winDAF[syn==4 & polym==TRUE]))
nsDAF <- as.vector(na.omit(winDAF[syn==0 & polym==TRUE]))

Ps.less5 <- sum(synDAF<0.05) # P0 MAF less than 5%
Ps.more5 <- sum(synDAF>0.05) # P0 MAF more than 5% 
Pns.less5 <- sum(nsDAF<0.05) # Pi MAF less than 5% 
Pns.more5 <- sum(nsDAF>0.05) # Pi MAF more than 5%

Pns.neutral.less5 <- Pns*(Ps.less5/Ps) # Proportion of neutral within the MAF < 5% class
Pns.neutral <- Pns.neutral.less5 + Pns.more5 # For alpha
Pns.weak <- Pns.less5 - Pns.neutral.less5

alpha.cor <- 1-(Pns.neutral/Ps)*(Ds/Dns)
contingency <- matrix(c(Pns.neutral,Ps,Dns,Ds),c(2,2))
test <- if(!is.na(sum(contingency))){
  if(sum(contingency)>0){fisher.test(contingency)$p.value}
} else {NA}

# To fully implement extended MKT, we need ms/mns, that is, the number of sites
# of each class. Doing that would require modifying the 'set.synnonsyn2' function
# to obtain all codons and check fold in every position (0,1,2)

## ADD NEW ROW ##
newrow <- c(S,Pi(k,m,n),DAF,divsites,D,K,unknowns,alpha.cor,test,Pns.neutral,Ps,Dns,Ds)
print(newrow)
tabsum[window,] <- newrow
}
return(tabsum)  
}
