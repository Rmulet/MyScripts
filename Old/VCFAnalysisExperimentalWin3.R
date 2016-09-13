#!/usr/bin/Rscript

#setwd("~/Documents/2_GenomicsData/TestPopGenome")
#filename <- "merge.vcf.gz"; ini <- 31768081-1; end <- 31899603; chrom <- "22"; wsize <- 10000; MK <- TRUE; ref.chr="chr22.fa"

# Region Analysis v0.8 - Imports the VCF file containing human and chimpanzee data and calculates
# metrics of polymorphism and divergence. The output is sent to a MySQL database.

# WARNING: The program is mean to be called from the terminal (VCFmerger.sh)
# UPDATE: 1) Tests of neutrality added; 2) MKT removed; 3) SFS is calculated from the derived allele.
# UPDATE2: Monomorphic positions in humans are not considered for SFS calculation. Mean is replaced with histogram.
# UPDATE3: Windows are converted into 0-based (BED format). Instead of a single column with the range, two columns (start,end) 
# are created. Generalization to multiple chromosomes is available.
# UPDATE4: DAF instead of SFS. Dots to underscores for MySQL. DB name changeable.
# UPDATE5: MKT re-added, optional. It requires a GFF and a FASTA file.
# UPDATE6: 1) Improved DAF with 1000 AA 2) Improved MKT with 4fould and extended MKT.

library(parallel)
no.cores <- detectCores() -1
c1 <- makeCluster(no.cores)

LocationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script
{
  this.file = NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
  }
  
  if (!is.null(this.file)) return(dirname(this.file))
  
  # But it may also be called from the command line
  cmd.args <- commandArgs(trailingOnly = FALSE)
  cmd.args.trailing <- commandArgs(trailingOnly = TRUE)
  cmd.args <- cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
  res <- gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
  
  # If multiple --file arguments are given, R uses the last one
  res = tail(res[res != ""], 1)
  if (0 < length(res)) return(dirname(res))
  
  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}
current.dir <- LocationOfThisScript()

#############################
## IMPORT AND PREPARE DATA ##
#############################

if (suppressMessages(!require("PopGenome"))) { # Can be removed
  print ("The 'PopGenome' package is missing and will be installed")
  install.packages("PopGenome")
  library("PopGenome")
}

args <- commandArgs(trailingOnly = TRUE) # Import arguments from command line
filename <- args[1] # Name of the file specified after the script
chrom <- args[2] # Chromosome number
ini <- as.numeric(args[3])-1; end <- args[4] # Window range. -1 from ini to include position 1.
wsize <- as.numeric(args[5]) # Window size
db <- args[6] # Name of the database [Genomics]
MK <- args[7] # Calculate MKT

# (1) We need to keep the .tbi file in the same folder as the vcf.gz (which must be compressed)
# (2) The chromosome identifier in the GFF has to be identical to the identifier used in the VCF file
# Since VCFs from the 1000GP use numeric-only identifiers (e.g. '22'), a conversion is required:
# In bash: cat Chr22.gff | sed 's/^chr22/22'. Also, tid must be identical to those two.
# (3) Numcols indicates the number of SNPs read into the RAM at once. For a sample of 1000 individuals,
# 10,000 are recommended on a 4 GB RAM computer.

init <- Sys.time()

if (MK==TRUE) {
  region <- readVCF(filename,numcols=5000,tid=chrom,from=ini,to=end,include.unknown=TRUE,gffpath=sprintf("chr%s.gff",chrom))
  # source(paste(current.dir,"set.synnonsyn_alt.R",sep="/"))
  source(paste("/home/roger/Documents/Scripts","set.synnonsyn_alt.R",sep="/")) # Only testing!!!
  # We use 'set.synnonsyn2' to identify four-fold degenerate
  region <- set.synnonsyn2(region,ref.chr=sprintf("chr%s.fa",chrom),save.codons=TRUE) 
  cat("MKT will be calculated \n")
} else {
  region <- readVCF(filename,numcols=5000,tid=chrom,from=ini,to=end,include.unknown=TRUE)
}

# Only runs on Linux, mcapply (used in parallelpackage) does not work on Windows
# region <- readVCF(filename,numcols=5000,parallel=TRUE,tid=chrom,from=ini,to=end,include.unknown=TRUE)

## DEFINE POPULATION/OUTGROUP ##

# The outgroup can be defined as a different population so that calculations of diversity and neutrality 
# are performed separately. 

individuals <- get.individuals(region)[[1]]
humans <- individuals[which(individuals != c("Chimp","Chimp.2"))]
chimp <- individuals[which(individuals == c("Chimp","Chimp.2"))]
region <- set.populations(region,list(humans,chimp))
region <- set.outgroup(region,chimp)

# The two steps above must be done BEFORE splitting data.

## SLIDING WINDOWS ##

# Sliding windows of size = wsize (1000 by default). Type=2 for genomic regions (type=1 for SNPs).

slide <- sliding.window.transform(region,width=wsize,jump=wsize,type=2) # Create windows in the object slide.
windows <- unlist(strsplit(slide@region.names,":")) # Retrieves window labels without ":"
windows <- t(sapply(strsplit(windows,"-"),as.numeric))
colnames(windows) <- c("start","end")
nwin <- nrow(windows)
# winsites <- slide@region.data@biallelic.sites # Positions in each window
allwindex <- slide@SLIDE.POS

################
## STATISTICS ##
################

## NEUTRALITY TESTS ##

# This calculates neutrality stats (Watterson's Theta, Fu and Li, Tajima's D) and the
# number of segregating sites (S). However, S excludes missing information in the outgroup.
slide <- neutrality.stats(slide)

## SITE FREQUENCY SPECTRA (IMPROVED) ##

n <- length(humans) # N = samples. Diploid, not divided by 2
# bial <- get.biallelic.matrix(region,1) # Testing
bialhuman <- get.biallelic.matrix(region,1) [1:n,,drop=F] #  Biallelic matrix, remove outgroup
freqs <- colSums(bialhuman == 1)/nrow(bialhuman)
freqs.df <- data.frame(POS=names(freqs),MAF=unname(freqs)) # Position and DAF of biallelic variants

# IDENTIFY ANCESTRAL ALLELES (AA) FOR SNPS IN THE REGION
library(stringr)
# temp <- system(sprintf('gunzip -c %s | cut -f2,4,8',filename),intern=TRUE) # Extract all variants from the merge.vcf.gz file  
# gpimport <- read.table(textConnection(temp),sep="\t",stringsAsFactors = FALSE,header=TRUE)
gpimport <- read.table("gpimport.txt",sep="\t",stringsAsFactors = FALSE,header=TRUE)
gpances <- toupper(sapply(gpimport[,3],function(x){str_match(x,"AA=([:alpha:])\\|")[2]},USE.NAMES=FALSE)) # Capture the ancestral allele
aarefs <- na.omit(cbind(gpimport[,1:2],ANC=gpances)) # Na.omit removes positions without known AA and divergent sites

# CALCULATE DERIVED ALLELE FREQUENCY (DAF)
mafan <- merge(freqs.df,aarefs) # Intersection between SNPs in the biallelic matrix and the AA in the VCF file
mafan$DAF <- NA # DAF is equal to MAF if derived = minor, but 1 - MAF if otherwise
mafan[mafan$REF != mafan$ANC,]$DAF <- 1- mafan[mafan$REF != mafan$ANC,]$MAF 
mafan[mafan$REF == mafan$ANC,]$DAF <- mafan[mafan$REF == mafan$ANC,]$MAF 
mafan <- merge(freqs.df[,1,drop=F],mafan,all.x=TRUE) # Retain originals
# Frequencies do not exactly match those of the VCF file! Why? Rounding perhaps
# Some sites do not have AA. In some cases because they are divergent, but others are polymorphic:
# freqs.df[which(is.na(mafan$MAF))[!(which(is.na(mafan$MAF)) %in% which(freqs.df$MAF == 0)),]

## CUSTOM FUNCTION FOR NUCLEOTIDE DIVERSITY AND DIVERGENCE ##
measures <- function(object) {
  n <- length(humans) # N = samples. Diploid, not divided by 2
  total <- length(individuals) # Total number of samples (including outgroup)
  
  ## POLYMORPHISM ## 
  # PopGenome has its own 'neutrality' and 'diversity' stats modules, but diversity takes a 
  # very long time (55 seconds including Pi for a region of 100 kb). 
  Pi <- function(k,m,n) { # Window number
    comb <- choose(n,2) # Binomial coefficient = combination without repetition
    pi <- k/(comb*m)
    return(round(pi,7))
  }
  # TABLE CONTAINING THE DATA:
  tabsum <- as.data.frame(matrix(numeric(nwin*13),ncol=13,nrow=nwin))
  colnames(tabsum) <- c("S","Pi","DAF","Divsites","D","K","Unknown","Alpha","Fisher","Pns","Ps","Dns","Ds")

  # LOOP TO READ THE WINDOWS IN SLIDE VARIABLE:
  for (window in 1:nrow(windows)) {
    # DETERMINE SEGREGATING SITES (EXCLUDING OUTGROUP)
    bial <- get.biallelic.matrix(slide,window) # Biallelic matrix of the window
    if (is.null(bial)||dim(bial)[2]==0) { # When no variants are detected
      newrow <- rep(0,13) # Empty rows
      tabsum[window,] <- newrow
      next
    }
    bialhuman <- bial[1:n,,drop=F] # Remove outgroup (drop = F to keep 1 dimension)
    misshuman <- colSums(is.na(bialhuman))>0 # Sites missing in humans (e.g. structural variants)
    polym <- apply(bialhuman,2,sum)>0 & !misshuman # Sites polymorphic in humans w/o missing
    bialhuman <- bialhuman[,polym,drop=FALSE] # Keep only polymorphic sites
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
    ## SITE FREQUENCY SPECTRUM (SFS) DERIVED ALLELE FREQUENCY (DAF) ##
    # The 'winsites' variable contains the coordinates of the biallelic variants in each window. 
    # By using this information, we can subset the DAF of the sites in this window contained in 'mafan'.
    windex <- allwindex[[window]]
    winDAF <- as.vector(na.omit(mafan$DAF[windex]))
    DAF <- hist(winDAF,seq(0.0,1,0.05),plot=F)$counts
    DAF <- paste(DAF,collapse = ";")

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

    ## MKT CALCULATION ##
    
    # STANDARD: Same result as MKT function if done with all syn.

    syn <- slide@region.data@synonymous[[window]] # Contains syn and non-syn positions
    Ps <- sum(polym[syn==4],na.rm=TRUE) # For all syn, >=1
    Pns <- sum(polym[syn==0],na.rm=TRUE)
    Ds <- sum(divtotal[syn==4],na.rm=TRUE)
    Dns <- sum(divtotal[syn==0],na.rm=TRUE)
    NI <- (Pns/Ps)/(Dns/Ds)
    alpha <- 1-NI
    DoS <- Dns/(Dns+Ds)-Pns/(Pns+Ps) # Direction of selection
    
    # ADAPTED:
    
    # 1) Use four-fold degenerated sites 2) Remove sites according to DAF
    
    # Using DAF results in a loss of information, as it is not available for some
    # variants. For MAF, simply retrieve from freqs.df.
    
    #winDAF <- mafan$DAF     #Testing for all region, REMOVE
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

## INTEGRATION OF NEUTRALITY, DIVERSITY AND DIVERGENCE METRICS ## 
 
regiondata <- measures(slide)
S2 <- slide@n.segregating.sites[,1] # Segretaging sites excluding unknowns
Tajima_D <- round(slide@Tajima.D[,1]/wsize,7) # 
FuLi_F <- round(slide@Fu.Li.F[,1]/wsize,7)
theta <- round(slide@theta_Watterson[,1]/wsize,7)
if (exists("S2")) {regiondata <- cbind(regiondata[,1:2],theta,S2,Tajima_D,FuLi_F,regiondata[,3:ncol(regiondata)])
  } else { regiondata <- cbind(regiondata[,1:2],theta=NA,S2=0,Tajima_D=0,FuLi_F=0,regiondata[,3:ncol(regiondata)]) }

print(Sys.time()-init)

######################
## DATA EXPORTATION ##
######################

# We add the chromosome info and the coordinates. For convenience in the next steps:
# a) Coordinates are converted to BED format: 0-based, subtract 1 from start column.
# b) Chromosome is expressed in "chrNN" format.
windows[,1] <- windows[,1]-1
export <- cbind(chr=rep(paste(c("chr",chrom),collapse=""),NROW(windows)),windows,regiondata)

suppressMessages(library(DBI))
suppressMessages(library(RMySQL))

con <- dbConnect(RMySQL::MySQL(),
                 user="roger", password="RM333",
                 dbname="PEGH", host="158.109.215.40")

first <- !dbExistsTable(con,db)

dbWriteTable(con,value=export,name=db,row.names=F,append=T)

#if (first == TRUE) # Remove if we want to concatenate various chromosomes
#  dbSendQuery(con,sprintf("ALTER TABLE %s CHANGE COLUMN start start VARCHAR(30);",db))
#  dbSendQuery(con,sprintf("ALTER TABLE %s ADD PRIMARY KEY (start);",db))
# }

on.exit(dbDisconnect(con))

print(Sys.time() - start.time)