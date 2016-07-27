#!/usr/bin/Rscript
# filename <- "merge.vcf.gz"; ini <- 31768006; end <- 31899628; wsize <- 200

# Region Analysis v0.6 - Imports the VCF file containing human and chimpanzee data and calculates
# metrics of polymorphism and divergence. The output is sent to a MySQL database.

start.time <- Sys.time()

#############################
## IMPORT AND PREPARE DATA ##
#############################

suppressMessages(require("PopGenome",quietly=TRUE))

args <- commandArgs() # Import arguments from command line
filename <- args[6] # Name of the file specified after the script
ini <- args[7]; end <- args[8] # Window range
wsize <- as.numeric(args[9]) # Window size
print(typeof(wsize))

# (1) The chromosome identifier in the GFF has to be identical to the identifier used in the VCF file
# Since VCFs from the 1000GP use numeric-only identifiers (e.g. '22'), a conversion is required:
# In bash: cat Chr22.gff | sed 's/^chr22/22/g'. Also, tid must be identical to those two. 
# (2) We need to keep the .tbi file in the same folder as the vcf.gz (which must be compressed)
region <- readVCF(filename,numcols=5000,tid="22",from=ini,to=end,gffpath="chr22.gff",include.unknown=TRUE)

## DEFINE POPULATION/OUTGROUP ##

# The outgroup can be defined as a different population so that calculations of diversity and neutrality 
# are performed separately. 

individuals <- get.individuals(region)[[1]]
humans <- individuals[which(individuals != c("Chimp","Chimp.2"))]
chimp <- individuals[which(individuals == c("Chimp","Chimp.2"))]
region <- set.populations(region,list(humans,chimp))

## VERIFY SYNONYMOUS/NON-SYNONYMOUS ##

region <- set.synnonsyn(region,ref.chr="chr22.fa") 

# The two steps above must be done BEFORE splitting data.

## SLIDING WINDOWS ##

# Sliding windows of size = 200. Type=2 for genomic regions (type=1 for SNPs)

slide <- sliding.window.transform(region,width=wsize,jump=wsize,type=2) # Create windows in the object slide.
windows <- slide@region.names # Window labels
windows <- trimws(unlist(strsplit(windows,":")))

################
## STATISTICS ##
################

# PopGenome has its own 'neutrality' and 'diversity' stats modules, but diversity takes a 
# very long time (55 seconds including Pi for a region of 100 kb)
measures <- function(object) {
  n <- length(humans) # N = samples. Diploid, not divided by 2
  total <- length(individuals) # Total number of samples (including outgroup)
  
  ## POLYMORPHISM ## 
  Theta <- function(S,m,n) { # S=Segregating sites; m=total sites; n=number of samples
    summ <- 0
    for (i in 1:(n-1)) {
      summ <- summ + 1/i
    }
    theta <- (S/m)/summ
    return(round(theta,7))
  }
  
  Pi <- function(k,m,n) { # Window number
    comb <- choose(n,2) # Binomial coefficient = combination without repetition
    pi <- k/(comb*m)
    return(round(pi,7))
  }
  
  # TABLE CONTAINING THE DATA
  tabsum <- data.frame(S=numeric(0),Theta=numeric(0),Pi=numeric(0),Divsites=numeric(0),D=numeric(0),K=numeric(0),Unknown=numeric(0)) 
  
  # LOOP TO READ THE WINDOWS IN SLIDE VARIABLE:
  for (window in 1:length(windows)) {
    # DETERMINE SEGREGATING SITES (EXCLUDING OUTGROUP)
    bial <- get.biallelic.matrix(slide,window) # Biallelic matrix
    if (is.null(bial)||dim(bial)[2]==0) { # When no variants are detected
      newrow <- c(0,0,0,0,0,0,0)
      tabsum[nrow(tabsum)+1,] <- newrow
      next
    }
    bialhuman <- bial[1:n,,drop=F] # Remove outgroup (drop = F to keep 1 dimension)
    mono <- apply(bialhuman,2,sum)>0 # Sites monomorphic in humans: all are REF alleles
    bialhuman <- bialhuman[,mono,drop=F] # Remove sites monomorphic in humans
    misshuman <- colSums(is.na(bialhuman)) # Sites missing in humans
    bialhuman <- bialhuman[,misshuman == 0,drop=FALSE] # Remove missing (not expected)
    # DETERMINE M (EXCLUDING MISSING AND POLYALLELIC)
    win <- as.numeric(strsplit(windows[window],"-")[[1]][1])
    # Total number of sites: 1) Remove NA in humans 2) Remove polyallelic sites
    if (!length(region@region.data@polyallelic.sites) == 0) { # Make sure list exists
      polyal <- region@region.data@polyallelic.sites[[1]]
      polysites <- sum(!is.na(match(polyal,win:(win+wsize))))
    }
    else {
      polysites <- 0
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
    # COUNT THE NUMBER OF UNKNOWNS IN THE OUTGROUP:
    unknowns <- sum(is.nan(bial[total,])) # Displayed as NaN
    
    ## DIVERGENCE ##
    
    # COUNT THE NUMBER OF DIVERGENCE SITES:
    divtotal <- bial[total,] == 1 # 1/1 in the outgroup (including polymorphic)
    divsites <- sum (divtotal & !mono,na.rm=T) # 1/1 in the outgroup (excluding polymorphic)
    # CALCULATE D AND K:
    mout <- m - unknowns
    D <- round(divsites/mout,7) # Observed divergence: Proportion of sites with divergent nucleotides
    K <- round(-3/4*log(1-4/3*D),7) # Real divergence: Jukes and Cantor model
  
    ## ADD NEW ROW ##
    newrow <- c(S,Theta(S,m,n),Pi(k,m,n),divsites,D,K,unknowns)
    tabsum[nrow(tabsum)+1,] <- newrow
  }
  return(tabsum)  
}

regiondata <- measures(slide)

## SITE FREQUENCY SPECTRA ##

# SFS is the default calculation in the module detail.stats. The results are stored in
# the slot slide@region.stats@minor.allele.freqs
slide <- detail.stats(slide)

# To calculate the mean SFS for each window:

SFSmean <- sapply(slide@region.stats@minor.allele.freqs, function(x){
  if(length(x)==0){return(0)}
  return(round(mean(x[1,], na.rm=TRUE),7))
  })

regiondata <- cbind(regiondata[,1:3],SFSmean,regiondata[,4:7]) # Reorder the fields

######################
## DATA EXPORTATION ##
######################

# We add the chromosome info and the windows coordinates:
export <- cbind(Chr=rep(22,length(windows)),Window=windows,regiondata)

suppressMessages(require(DBI))
suppressMessages(require(RMySQL))

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