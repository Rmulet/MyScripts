#!/usr/bin/Rscript
#setwd("~/Documents/2_GenomicsData/TestPopGenome")
#filename <- "merge.vcf.gz"; ini <- 31768006; end <- 31899628; wsize <- 1000

# Region Analysis v0.7 - Imports the VCF file containing human and chimpanzee data and calculates
# metrics of polymorphism and divergence. The output is sent to a MySQL database.
# WARNING: The program is mean to be called from the terminal (VCFmerger.sh)
# UPDATE: 1) Tests of neutrality added; 2) MKT removed; 3) SFS is calculated from the derived allele.
# UPDATE2: Monomorphic positions in humans are not considered for SFS calculation. Mean is replaced with histogram.
# UPDATE3: Windows are converted into 0-based (BED format) to facilitate the following analyses. Instead of a single
# column with the range, two columns (start,end) are created.

start.time <- Sys.time()

#############################
## IMPORT AND PREPARE DATA ##
#############################

suppressMessages(library("PopGenome",quietly=TRUE))

args <- commandArgs() # Import arguments from command line
filename <- args[6] # Name of the file specified after the script
ini <- args[7]; end <- args[8] # Window range
wsize <- as.numeric(args[9]) # Window size

# We need to keep the .tbi file in the same folder as the vcf.gz (which must be compressed)
# Optionally, we can provide a GFF file with annotations, but it's not necessary here.
# Numcols indicates the number of SNPs read into the RAM at once. For a sample of 1000 individuals,
# 10,000 are recommended on a 4 GB RAM computer.
region <- readVCF(filename,numcols=5000,tid="22",from=ini,to=end,include.unknown=TRUE)

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
windows <- trimws(do.call(rbind,strsplit(windows,"-"))) # Separates in start and end columns
windows <- apply(windows,2,as.numeric) # Makes the windows variable numeric

################
## STATISTICS ##
################

## NEUTRALITY TESTS ##

# This calculates neutrality stats (Watterson's Theta, Fu and Li, Tajima's D) and the
# number of segregating sites (S). However, S excludes missing information in the outgroup.
slide <- neutrality.stats(slide)

## SITE FREQUENCY SPECTRA ##

# SFS is the default calculation in the module detail.stats. The results are stored in
# the slot slide@region.stats@minor.allele.freqs
# slide <- detail.stats(slide)
# Since it does not work with an outgroup, we can try the following:
region <- detail.stats(region)
totalMAF <- region@region.stats@minor.allele.freqs
windex <- slide@SLIDE.POS # Index of each variant grouped by windows
MAFsites <- as.numeric(colnames(totalMAF[[1]])) # Only when outgroup is defined. If not, totalMAF does not have names.

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
  # TABLE CONTAINING THE DATA
  tabsum <- data.frame(S=numeric(0),Pi=numeric(0),SFS=numeric(0),Divsites=numeric(0),D=numeric(0),K=numeric(0),Unknown=numeric(0)) 

  # LOOP TO READ THE WINDOWS IN SLIDE VARIABLE:
  for (window in 1:nrow(windows)) {
    # DETERMINE SEGREGATING SITES (EXCLUDING OUTGROUP)
    bial <- get.biallelic.matrix(slide,window) # Biallelic matrix
    if (is.null(bial)||dim(bial)[2]==0) { # When no variants are detected
      newrow <- c(0,0,0,0,0,0,0)
      tabsum[nrow(tabsum)+1,] <- newrow
      next
    }
    bialhuman <- bial[1:n,,drop=F] # Remove outgroup (drop = F to keep 1 dimension)
    polym <- apply(bialhuman,2,sum)>0 # Sites polymorphic in humans: non REF (0) alleles
    bialhuman <- bialhuman[,polym,drop=F] # Remove sites monomorphic in humans
    misshuman <- colSums(is.na(bialhuman)) # Sites missing in humans
    bialhuman <- bialhuman[,misshuman == 0,drop=FALSE] # Remove missing (not expected)
    # DETERMINE M (EXCLUDING MISSING AND POLYALLELIC)
    winstart <- windows[window,1] # Select the start position in that window
    # Total number of sites: 1) Remove NA in humans 2) Remove polyallelic sites
    if (!length(region@region.data@polyallelic.sites) == 0) { # Make sure list exists
      polyal <- region@region.data@polyallelic.sites[[1]] # Positions of all polyalleles
      polysites <- sum(!is.na(match(polyal,winstart:(winstart+wsize)))) # N polyallelic sites in that window
    }
    else {
      polysites <- 0 # If not available, we asumme 0
    }
    m <- wsize-sum(misshuman,na.rm=T)-polysites
    # DETERMINE S AND K (WHEN VARIANTS ARE AVAILABLE)
    if (is.null(bialhuman)||dim(bialhuman)[2]==0) {
      S <- 0
      k <- 0
    }
    else {
      S <- length(bialhuman[1,]) # Number of variants
      freqs <- apply(bialhuman,2,table)
      k <- sum(freqs[1,]*freqs[2,]) # Note that k and K are different!
    }
    # CALCULATE SITE FREQUENCY SPECTRUM (SFS):
    # When there is an outgroup, missing positions are eliminated in MAF, but not in the biallelic matrix:
    # the indices of the windows do not match those in the totalMAF variable. Instead, we have to 
    # use the names of the variable (MAFsites). We filter to exclude monomorphic sites.
    if (length(MAFsites) != 0) { # Outgroup: indices of 'totalMAF' =/= those of windows.
      winsites <- which(MAFsites %in% windex[[window]][polym]) # Window sites w/o unknowns and polymorphic in humans
      }
    else { # No outgroup: indices of 'totalMAF' match those of windows in 'windex'- 
      winsites <- windex[[window]][polym] # Sites in that window w/o unknowns and polymorphic in humans
    }
    winMAF <- totalMAF[[1]][1,winsites] # Frequencies in window (removing unknowns)
    SFS <- hist(winMAF,seq(0.0,1,0.1),plot=F)$counts
    SFS <- paste(SFS,collapse = ";")
    # COUNT THE NUMBER OF UNKNOWNS IN THE OUTGROUP:
    unknowns <- sum(is.nan(bial[total,])) # Displayed as NaN in the biallelic matrix
    
    ## DIVERGENCE ##
    
    # COUNT THE NUMBER OF DIVERGENCE SITES:
    divtotal <- bial[total,] == 1 # 1/1 in the outgroup (including polymorphic)
    divsites <- sum (divtotal & !polym,na.rm=T) # 1/1 in the outgroup (excluding polymorphic)
    # CALCULATE D AND K:
    mout <- m - unknowns
    D <- round(divsites/mout,7) # Observed divergence: Proportion of sites with divergent nucleotides
    K <- round(-3/4*log(1-4/3*D),7) # Real divergence: Jukes and Cantor model
  
    ## ADD NEW ROW ##
    newrow <- c(S,Pi(k,m,n),SFS,divsites,D,K,unknowns)
    tabsum[nrow(tabsum)+1,] <- newrow
  }
  return(tabsum)  
}

## INTEGRATION OF NEUTRALITY, DIVERSITY AND DIVERGENCE METRICS ## 
 
regiondata <- measures(slide)
S2 <- slide@n.segregating.sites[,1]/wsize # Segretaging sites excluding unknowns
Tajima.D <- slide@Tajima.D[,1]/wsize
FuLi.F <- slide@Fu.Li.F[,1]/wsize
theta <- slide@theta_Watterson[,1]/wsize
regiondata <- cbind(regiondata[,1:2],theta,S2,Tajima.D,FuLi.F,regiondata[,3:7])

######################
## DATA EXPORTATION ##
######################

# We add the chromosome info and the windows coordinates. Convert to BED format, which is 0-based. 
export <- cbind(Chr=rep("chr22",nrow(windows)),Window=windows[,1]-1,regiondata)

suppressMessages(library(DBI))
suppressMessages(library(RMySQL))

con <- dbConnect(RMySQL::MySQL(),
                 user="root", password="RM333",
                 dbname="PEGH", host="localhost")

first <- !dbExistsTable(con,"Genomics")

dbWriteTable(con,value=export,name="Genomics",row.names=F,append=T)

if (first == TRUE) {
  dbSendQuery(con,"ALTER TABLE Genomics CHANGE COLUMN Window Window VARCHAR(30);")
  dbSendQuery(con,"ALTER TABLE Genomics ADD PRIMARY KEY (Window);")
  }

on.exit(dbDisconnect(con))

print(Sys.time() - start.time)