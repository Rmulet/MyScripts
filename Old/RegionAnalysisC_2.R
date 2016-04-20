require("PopGenome")
setwd("~/Documents/2_GenomicsData/TestPopGenome")

start.time <- Sys.time()

## IMPORT DATA ##

# (1) The chromosome identifier in the GFF has to be identical to the identifier used in the VCF file
# Since VCFs from the 1000GP use numeric-only identifiers (e.g. '22'), a conversion is required:
# In bash: cat Chr22.gff | sed 's/^chr22/22'. Also, tid must be identical to those two. 
# (2) We need to keep the .tbi file in the same folder as the vcf.gz (which must be compressed)
region <- readVCF("merge.vcf.gz",numcols=5000,tid="22",from=31768006,to=31799567,gffpath="chr22.gff",include.unknown=TRUE)

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

slide <- sliding.window.transform(region,width=200,jump=200,type=2) # Create windows in the object slide.
windows <- slide@region.names # Window labels
windows <- trimws(unlist(strsplit(windows,":")))

################
## STATISTICS ##
################

## POLYMORPHISM ##

# PopGenome has its own 'neutrality' and 'diversity' stats modules, but this is faster.
polymorphism <- function(object) {
  n <- length(humans) # N = samples. Diploid, not divided by 2
  total <- length(individuals) # Total number of samples (including outgroup)
  
  Theta <- function(S,m,n) { # S=Segregating sites; m=total sites; n=number of samples
    summ <- 0
    for (i in 1:(n*2-1)) {
      summ <- summ + 1/i
    }
    theta <- (S/m)/summ
    return(theta)
  }
  
  Pi <- function(k,m,n) { # Window number
    comb <- choose(n,2) # Binomial coefficient = combination without repetition
    pi <- k/(comb*m)
    return(pi)
  }
  
  # TABLE CONTAINING THE DATA
  tabsum <- data.frame(S=numeric(0),Theta=numeric(0),Pi=numeric(0),Unknown=numeric(0),Divsites=numeric(0)) 
  
  # LOOP TO READ THE WINDOWS IN SLIDE VARIABLE:
  for (window in 1:length(windows)) {
  # DETERMINE SEGREGATING SITES (EXCLUDING OUTGROUP)
    bial <- get.biallelic.matrix(slide,window) # Biallelic matrix
    if (is.null(bial)||dim(bial)[2]==0) { # When no variants at all are detected
      newrow <- c(0,0,0,0,0)
      tabsum[nrow(tabsum)+1,] <- newrow
      rownames(tabsum)[window] <- windows[window]
      next
    }
    bialh <- bial[1:n,,drop=F] # Remove outgroup
    bialhuman <- bialh[,(apply(bialh,2,sum)>0),drop=F] # Remove SNPs monomorphic in humans
    bialhuman <- bialhuman[,colSums(is.na(bialhuman)) == 0,drop=FALSE] # Remove missing (not expected)
    # DETERMINE M (EXCLUDING MISSING AND POLYALLELIC)
    win <- as.numeric(strsplit(windows[window],"-")[[1]][1])
    # Total number of sites: 1) Remove NA in humans 2) Remove polyallelic sites
    if (!length(region@region.data@polyallelic.sites) == 0) { # Make sure list exists
      polyal <- region@region.data@polyallelic.sites[[1]]
      polysites <- sum(!is.na(match(polyal,win:(win+200))))
    }
    else {
      polysites <- 0
    }
      m <- slide@n.sites[window]-sum(colSums(is.na(bialhuman)))-polysites
    # DETERMINE S AND K (WHEN VARIANTS ARE AVAILABLE)
    if (is.null(bialhuman)||dim(bialhuman)[2]==0) {
      S <- 0
      k <- 0
    }
    else {
      S <- length(bialhuman[1,]) # Number of variants ()
      freqs <- apply(bialhuman,2,table)
      k <- sum(freqs[1,]*freqs[2,])
    }
    # COUNT THE NUMBER OF UNKNOWNS IN THE OUTGROUP:
    unknowns <- sum(is.nan(bial[total,])) # Displayed as NaN
    # COUNT THE NUMBER OF DIVERGENCE SITES:
    divsites <- sum(bial[total,] == 1,na.rm = T) # 1/1 in the outgroup (not excluding polymorphic)
    # ADDING A NEW ROW TO THE TABLE:
    newrow <- c(S,Theta(S,m,n),Pi(k,m,n),unknowns,divsites)
    tabsum[nrow(tabsum)+1,] <- newrow
    rownames(tabsum)[window] <- windows[window]
  }
  return(tabsum)  
}

regiondata <- polymorphism(slide)

# SFS is the default calculation in the module detail.stats. The results are stored in
# the slot slide@region.stats@minor.allele.freqs
slide <- detail.stats(slide)

# To calculate the mean SFS for each window:

SFSmean <- sapply(slide@region.stats@minor.allele.freqs, function(x){
if(length(x)==0){return(0)}
return(mean(x[1,], na.rm=TRUE))
})

regiondata <- cbind(regiondata,SFSmean)

## DIVERGENCE ##

# MKT is automatically conducted for every pair of populations. The outgroup must be defined as
# a population with set.populations (no set.outgroup necessary). Synonyms are a required!

slide <- sliding.window.transform(region,width=200,jump=200,type=2) # Create windows in the object slide.
slide <- MKT(slide)
colnames(slide@MKT)[c(5,7)] <- c("NI","p-val")
regiondata <- cbind(regiondata,slide@MKT[,1:6])

# If only for genes:

#split <- splitting.data(region,subsites="gene")
#split <- MKT(split); split@MKT

# Or the entire region: region <- MKT(region);region@MKT

#######################
## DATA EXPOORTATION ##
#######################

export <- regiondata[,c(1,2,3,6,5,7,8,9,10,11,4)]
Sys.time() - start.time

require(DBI)
require("RMySQL")
con <- dbConnect(RMySQL::MySQL(),
                 user="root", password="RM333",
                 dbname="PEGH", host="localhost")
dbWriteTable(con,value=export,name="Genomics",overwrite=T,row.names=F,col.names=F)

on.exit(dbDisconnect(con))