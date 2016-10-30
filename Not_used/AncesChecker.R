#!/usr/bin/Rscript
#setwd("/home/roger/Documents/4_Analysis/Genomics")

start.time <- Sys.time()

#############################
## IMPORT AND PREPARE DATA ##
#############################

if (suppressMessages(!require("PopGenome"))) { # Can be removed
  print ("The 'PopGenome' package is missing and will be installed")
  install.packages("PopGenome")
  library("PopGenome")
}

filename <- "ancesinspect.vcf.gz"; ini <- 20000000; end <- 20010000; wsize <- 1000

args <- commandArgs(trailingOnly = TRUE) # Import arguments from command line
filename <- args[1] # Name of the file specified after the script
chrom <- args[2] # Chromosome number
ini <- as.numeric(args[3])-1; end <- args[4] # Window range. -1 from ini to include position 1.
wsize <- as.numeric(args[5]) # Window size

selregion <- function(filename="inspect.vcf.gz",ini=20000000,end=20100000) {
  setwd("/home/roger/Documents/4_Analysis/Genomics")
  system(sprintf("bcftools merge -Oz --missing-to-ref -o %s -r 22:%d-%d chr22_filtered.vcf.gz chr22_aln.vcf.gz",filename,ini,end))
  system(sprintf("tabix %s",filename))
}
  
selregion(filename)
region <- readVCF(filename,numcols=5000,tid="22",from=ini,to=end,include.unknown=TRUE)

## DEFINE POPULATION/OUTGROUP ##

# The outgroup can be defined as a different population so that calculations of diversity and neutrality 
# are performed separately. 

individuals <- get.individuals(region)[[1]]
humans <- individuals[which(individuals != c("Chimp","Chimp.2"))]
chimp <- individuals[which(individuals == c("Chimp","Chimp.2"))]
region <- set.populations(region,list(humans,chimp))
region <- set.outgroup(region,chimp)

############################
## SITE FREQUENCY SPECTRA ##
############################

## CALCULATE FREQUENCIES ##

# SFS is the default calculation in the module detail.stats. The results are stored in
# the slot slide@region.stats@minor.allele.freqs
region <- detail.stats(region)
totalMAF <- region@region.stats@minor.allele.freqs # MAF of all the loci in the window
MAFh <- totalMAF[[1]][1,] # MAF of the human population

## SNPS WITH FREQ > 0.9

freqs = MAFh[MAFh > 0.9 & MAFh < 1] # SNPs with MAF superior to 0.9 but inferior to 1
mafsites = as.numeric(names(freqs)) # Indices of SNPs with MAF 0.9-1
bial <- get.biallelic.matrix(region,1) # Biallelic matrix
bialfil = bial[,mafsites] # Biallelic matrix filtered for SNPs with MAF 0.9-1
ancsnps <- colnames(bialfil) # Positions of those SNPs
write(ancsnps,"SFStest.txt")

temp <- system('for string in `cat SFStest.txt`; do grep "$string" merge.vcf;done | cut -f1-8,2271',intern=TRUE)
ances <- read.table(textConnection(temp),sep="\t",stringsAsFactors = FALSE)
# do.call(rbind,strsplit(ances[,8],";"))[,c(1,2)]
library(stringr)
equals= 0
for (i in 1:nrow(ances)) {
  gpancestral <- str_match(ances[i,8],"AA=(.)\\|")[2]
  chimp <- ances[i,5] # Assumed ancestral = alternative allele (all chimps)
  print(c(gpancestral,chimp))
  if (gpancestral == chimp) {equals = equals+1}
}
print((1-equals/nrow(ances))*100) # Estimation of the proportion of wrongly assigned ancestral alleles
wrong <- nrow(ances)
right <- equals

## REST OF THE SNPS (freq 0-0.9)

freqs = MAFh[MAFh < 0.9] # SNPs with MAF superior to 0.9 but inferior to 1
mafsites2 = as.numeric(names(freqs)) # Indices of SNPs with MAF 0.9-1
bial <- get.biallelic.matrix(region,1) # Biallelic matrix
bialfil2 = bial[,mafsites2] # Biallelic matrix filtered for SNPs with MAF 0.9-1
ancsnps2 <- colnames(bialfil2) # Positions of those SNPs
write(ancsnps2,"SFStest2.txt")

temp <- system('for string in `cat SFStest2.txt`; do grep "$string" merge.vcf;done | cut -f1-8,2271',intern=TRUE)
ances2 <- read.table(textConnection(temp),sep="\t",stringsAsFactors = FALSE)
library(stringr)
equals2= 0
for (i in 1:nrow(ances2)) {
  gpancestral <- str_match(ances2[i,8],"AA=(.)\\|")[2]
  if (is.na(gpancestral)) {next}
  if (ances2[i,9] == "1/1") {
    chimp <- ances2[i,5] # The assumed ancestral is the alternative allele
  } else {chimp <- ances2[i,4]} # The assumed ancestral is the ref allele
    print(c(gpancestral,chimp))
  if (gpancestral == chimp) {
    equals2 = equals2+1
  }
}
print((1-equals2/nrow(ances2))*100) # Estimation of the proportion of wrongly assigned ancestral alleles
wrong2 <- nrow(ances2)-equals2
right2 <- equals2

## Chi-Square test to verify it:

companc <- data.frame(Less0.9=c(right,wrong),All=c(right2,wrong2),row.names=c("right","wrong"))
chisq.test(companc,Yates=FALSE)














# TotalMAF already takes into account the presence of an outgroup in the entire window
n <- length(humans) # N = samples. Diploid, not divided by 2
bial <- get.biallelic.matrix(region,1) # Biallelic matrix
bialhuman <- bial[1:n,,drop=F] # Remove outgroup (drop = F to keep 1 dimension)
polym <- apply(bialhuman,2,sum)>0 # Sites polymorphic in humans: non REF (0) alleles
polybial <- bial[,polym]

# Some alleles in polybial do not have frequencies in MAFsites.Why?

sum(apply(bialhuman,2,sum)>0) # 281 polymorphic sites
seq(1,387)[!seq(1,387) %in% as.numeric(names(MAFh))]

polyMAF <- totalMAF[[1]][1,][polym] # We remove monomorphic sites with freq = 1 (fixed)
DAF <- hist(polyMAF,seq(0.0,1,0.1),plot=F)$counts

## SLIDING WINDOWS ##

windex <- slide@SLIDE.POS # Index of each variant grouped by windows
MAFsites <- as.numeric(colnames(totalMAF[[1]])) # Only when outgroup is defined. If not, totalMAF does not have names.