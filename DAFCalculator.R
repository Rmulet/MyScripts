#!/usr/bin/Rscript

# setwd("/home/roger/Documents/Genomics/MyScripts")
# install.packages("PopGenome",repos="http://cran.r-project.org")
# ini=45886885; end=45920909; chrom <- "22"; filename <- "chr22.vcf.gz"

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
      DAFCalculator - A simple script to calculate DAF using 1000 GP data and PopGenome.

      Usage: DAFCalculator.R <VCF-file> <chr> <start> <end> [-o|--output] [-h]

      Positional arguments:
      <vcf-file> STR            A VCF file containing the variants of interest
      <chr> INT                 Chromosome of the region under study (e.g. 22)
      <start> INT               Start of the region under study (e.g. 1)
      <end> INT                 End of the region under study (e.g. 1500)

      Optional arguments:
      --interval INT            Interval of frequencies for the DAF histogram
      --output FILENAME         Name of the output [DEFAULT: INPUT_daf.txt]
      -h, --help                Print this text

      Example usage: Rscript DAFcal.R chr21.vcf.gz 21 100 20000
      ")
  q(save="no") # Exit without saving
}

## Retrieve positional arguments
filename <- args[1]

if(!file.exists(filename)) {
  cat("\nFile does not exist! Exiting\n")
  q(save="no")
}

chrom <- args[2] # Chromosome number
ini <- as.numeric(args[3])-1; end <- as.numeric(args[4])-1 # Although the mask file is supposed to be in BED format,
# comparison with the FASTA file reveals that is it actually 1-based, half-open. That is, the end is not
# included. START matches the position in the FASTA, but PopGenome tends to add 1; therefore, we remove it

## Retrieve optional arguments

output.name <- "DAFTable"; interval <- 0.05

if(length(args) > 5) {
  optional.args <- args[5:length(args)]
  if (length(index <- grep("--output",optional.args)) != 0){
    output.name <- optional.args[index+1]
  }
  if (length(index <- grep("--interval",optional.args)) != 0){
    interval <- as.numeric(optional.args[index+1])
  }
}
  
cat(sprintf("Settings: %s:%s-%s %s %s\n",chrom,ini,end,output.name,interval))

########################
## LOAD GENOME OBJECT ##
########################

cat("Loading Popgenome object...\n")
suppressMessages(library("PopGenome"))

# We attempt to read the VCF file. This is not possible when, for instance, all human positions are monomorphic and some
# chimpanzee ones are unknown, probably due to an alignment with multiple sites at once. Another possibility is the complete
# absence of variants in that region, both between humans and chimpanzees and among humans.

region <- tryCatch({readVCF(filename,numcols=10000,tid=chrom,from=ini,to=end,include.unknown=TRUE)},error = function(e) {
  message(e); write(sprintf("%s:%d-%d -- %s",chrom,ini,end,e),sprintf("error_chr%s_%s.log",chrom,pop),append=TRUE)
  return(NULL)})

# Verify that the region object is not null (failure to load), contains variants and has been loaded onto R. If it is, then we assume
# that there are no variants, and therefore S,D and all related metrics are 0 (except for alpha)

if (is.null(region)||is.logical(region)||region@n.biallelic.sites==0) { # If readVCF fails, region=FALSE. If no variants, region=NULL
  cat("This region does not contain any variants: DAF set to NA\n\n")
  export <- data.frame(chr=chrom,start=ini,end=end,DAF=NA)
  write.table(export,file=paste(output.name,".tab",sep=""),quote=FALSE,sep="\t",row.names=F,append=TRUE,
              col.names=!file.exists(paste(output.name,".tab",sep="")))  # Column names written if file does not exist
  quit()
}

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

cat("\nCalculating minor allele frequencies...")
n <- length(humans) # N = samples. Diploid, not divided by 2
bialhuman <- get.biallelic.matrix(region,1) [1:n,,drop=F] #  Whole biallelic matrix, remove outgroup
MAF <- colSums(bialhuman == 1)/nrow(bialhuman)
MAF.df <- data.frame(POS=names(MAF),MAF=unname(MAF)) # Position and MAF of biallelic variants
cat("done\n")

## IDENTIFY ANCESTRAL ALLELES (AA) FOR SNPS IN THE REGION

cat("Obtaining ancestral allele from 1000KG file...")
library(stringr)
temp <- system(sprintf('gunzip -c %s | cut -f2,4,8',filename),intern=TRUE) # Extract all variants from the merge.vcf.gz file
gpimport <- read.table(textConnection(temp),sep="\t",stringsAsFactors = FALSE,header=TRUE)
# gpimport <- read.table("gpimport.txt",sep="\t",stringsAsFactors = FALSE,header=TRUE)
gpances <- toupper(sapply(gpimport[,3],function(x){str_match(x,"AA=([:alpha:])\\|")[2]},USE.NAMES=FALSE)) # Capture the ancestral allele
aarefs <- na.omit(cbind(gpimport[,1:2],ANC=gpances)) # Na.omit removes positions without known AA and divergent sites
cat("done\n")

## CALCULATE DERIVED ALLELE FREQUENCY (DAF)

cat("Calculating DAF based on ancestral allele...")
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
  DAF <- hist(mafan$DAF,seq(0.0,1,interval),plot=F)$counts
  DAF <- paste(DAF,collapse = ";")
} else {DAF <- NA}

cat("done\n")

######################
## DATA EXPORTATION ##
######################

export <- data.frame(chr=chrom,start=ini,end=end,DAF=DAF)
write.table(export,file=paste(output.name,".tab",sep=""),quote=FALSE,sep="\t",row.names=F,append=TRUE,
            col.names=!file.exists(paste(output.name,".tab",sep="")))  # Column names written if file does not exist
cat(sprintf("Data saved to file %s\n\n",paste(output.name,".tab",sep="")))
