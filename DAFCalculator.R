#!/usr/bin/Rscript

# setwd("D:/")
# install.packages("PopGenome",repos="http://cran.r-project.org")

#########################
## GET AND CHECK INPUT ##
#########################

## Collect arguments
args <- commandArgs(trailingOnly = TRUE)

## Default setting when no arguments passed
if(length(args) < 4) {
  args <- c("--help")
}

## Help section
if("--help" %in% args || "-h" %in% args) {
  cat("
      DAFCal - A simple script to calculate DAF using 1000 GP data and PopGenome.
      
      Usage: DAFCal.R <VCF-file> <chr> <start> <end> [-o|--output] [-h]

      Positional arguments:
      <vcf-file> STR            A VCF file containing the variants of interest
      <chr> INT                 Chromosome of the region under study (e.g. 22)  
      <start> INT               Start of the region under study (e.g. 1) 
      <end> INT                 End of the region under study (e.g. 1500)

      Optional arguments:
      -o, --output FILENAME     Name of the output [DEFAULT: INPUT_daf.txt]
      -h, --help                Print this text

      Example usage: Rscript DAFcal.R chr21.vcf.gz 21 100 20000
      ")
  q(save="no") # Exit without saving
}

## Get and check VCF filename
filename <- args[1]

if(!file.exists(filename)) {
  cat("\nFile does not exist! Exiting\n")
  q(save="no")
}

chrom <- args[2] # Chromosome number
ini <- as.numeric(args[3])-1; end <- as.numeric(args[4])-1 # Although the mask file is supposed to be in BED format,
# comparison with the FASTA file reveals that is it actually 1-based, half-open. That is, the end is not 
# included. START matches the position in the FASTA, but PopGenome tends to add 1; therefore, we remove it
 
if(args >= 5 & ("--output" %in% args |"-o" %in% args)) {
  output <- args[6]
}

########################
## LOAD GENOME OBJECT ##
########################

library("PopGenome")

region <- readVCF("chr22_olga.vcf.gz", numcols=10000, tid=chrom, from=ini, to=end, include.unknown = TRUE)

## DEFINE POPULATION/OUTGROUP ##
individuals <- get.individuals(region)[[1]]
humans <- individuals[which(individuals != c("Chimp","Chimp.2"))]
chimp <- individuals[which(individuals == c("Chimp","Chimp.2"))]
region <- set.populations(region,list(humans,chimp))
region <- set.outgroup(region,c("Chimp"), diploid=TRUE)

###################
## CALCULATE DAF ##
###################

## OBTAIN MAF

n <- length(humans) # N = samples. Diploid, not divided by 2
bialhuman <- get.biallelic.matrix(region,1) [1:n,,drop=F] #  Whole biallelic matrix, remove outgroup
MAF <- colSums(bialhuman == 1)/nrow(bialhuman)
MAF.df <- data.frame(POS=names(MAF),MAF=unname(MAF)) # Position and MAF of biallelic variants

## IDENTIFY ANCESTRAL ALLELES (AA) FOR SNPS IN THE REGION
library(stringr)
temp <- system(sprintf('gunzip -c %s | cut -f2,4,8',filename),intern=TRUE) # Extract all variants from the merge.vcf.gz file
gpimport <- read.table(textConnection(temp),sep="\t",stringsAsFactors = FALSE,header=TRUE)
# gpimport <- read.table("gpimport.txt",sep="\t",stringsAsFactors = FALSE,header=TRUE)
gpances <- toupper(sapply(gpimport[,3],function(x){str_match(x,"AA=([:alpha:])\\|")[2]},USE.NAMES=FALSE)) # Capture the ancestral allele
aarefs <- na.omit(cbind(gpimport[,1:2],ANC=gpances)) # Na.omit removes positions without known AA and divergent sites

## CALCULATE DERIVED ALLELE FREQUENCY (DAF)
mafan <- merge(MAF.df,aarefs) # Intersection between SNPs in the biallelic matrix and the AA in the VCF file
if (dim(mafan)[1] > 0) { # Avoid windows where variants have no known AA
  mafan$DAF <- NA # DAF is equal to MAF if derived = minor, but 1 - MAF if otherwise
  mafan[mafan$REF != mafan$ANC,]$DAF <- 1- mafan[mafan$REF != mafan$ANC,]$MAF
  mafan[mafan$REF == mafan$ANC,]$DAF <- mafan[mafan$REF == mafan$ANC,]$MAF
}
mafan <- merge(MAF.df[,1,drop=F],mafan,all.x=TRUE) # Retain originals

# Frequencies do not exactly match those of the VCF file! Why? Rounding perhaps
# Some sites do not have AA. In some cases because they are divergent, but others are polymorphic:
# MAF.df[which(is.na(mafan$MAF))[!(which(is.na(mafan$MAF)) %in% which(MAF.df$MAF == 0)),]

if (is.numeric (mafan$DAF)) {
  DAF <- hist(mafan$DAF,seq(0.0,1,0.05),plot=F)$counts
  DAF <- paste(DAF,collapse = ";")
} else {DAF <- NA}

######################
## DATA EXPORTATION ##
######################

export <- c(chr=chrom,start=ini,end=end,DAF=DAF)
write.table(export,file=paste(db,".tab",sep=""),quote=FALSE,sep="\t",row.names=F,append=TRUE,
            col.names=!file.exists(paste(db,".tab",sep="")))  # Column names written if file does not exist
