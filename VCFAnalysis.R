#!/usr/bin/Rscript

#setwd("~/Documents/2_GenomicsData/Final")
# filename <- "merge.3.vcf.gz"
#filename <- "merge.vcf.gz"; ini <- 31768081-1; end <- 31899603-1; chrom <- "22"; wsize <- 10000; MK <- TRUE; ref.chr="chr22.fa"; window = 7

# Region Analysis v0.8 - Imports the VCF file containing human and chimpanzee data and calculates
# metrics of polymorphism and divergence. The output is sent to a MySQL database.

# WARNING: The program is mean to be called from the terminal (VCFmerger.sh)
# UPDATE: 1) Tests of neutrality added; 2) MKT removed; 3) SFS is calculated from the derived allele.
# UPDATE2: Monomorphic positions in humans are not considered for SFS calculation. Mean is replaced with histogram.
# UPDATE3: Windows are converted into 0-based (BED format). Instead of a single column with the range, two columns (start,end)
# are created. Generalization to multiple chromosomes is available.
# UPDATE4: DAF instead of SFS. Dots to underscores for MySQL. DB name changeable.
# UPDATE5: MKT re-added, optional. It requires a GFF and a FASTA file.
# UPDATE6: 1) Improved DAF with 1000 AA 2) Improved MKT with 4fold and extended MKT.
# UPDATE7: GFFtoFASTA has replaced 'set.synnonsyn_alt' as the method of choic for MKT.
# UPDATE8: Change in the coordinate system according to the mask in pseudo-BED format.

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
ini <- as.numeric(args[3])-1; end <- as.numeric(args[4])-1 # Although the mask file is supposed to be in BED format,
# comparison with the FASTA file reveals that is it actually 1-based, half-open. That is, the end is not 
# included. START matches the position in the FASTA, but PopGenome tends to add 1; therefore, we remove it
wsize <- as.numeric(args[5]) # Window size
db <- args[6] # Name of the database [Genomics]
MK <- args[7] # Calculate MKT

# (1) We need to keep the .tbi file in the same folder as the vcf.gz (which must be compressed)
# (2) The chromosome identifier in the GFF has to be identical to the identifier used in the VCF file
# Since VCFs from the 1000GP use numeric-only identifiers (e.g. '22'), a conversion is required:
# In bash: cat Chr22.gff | sed 's/^chr22/22'. Also, tid must be identical to those two.
# (3) Numcols indicates the number of SNPs read into the RAM at once. For a sample of 1000 individuals,
# 10,000 are recommended on a 4 GB RAM computer.

start.time <- Sys.time()

region <- readVCF(filename,numcols=5000,tid=chrom,from=ini,to=end,include.unknown=TRUE)
# Syn-nonsyn is not needed if we only use 0- and 4-fold. Therefore, GFF and FASTA can be skipped.

# Verify that the region contains variants and has been loaded onto R.
if (region@n.biallelic.sites==0)||is.logical(region)) { # If readVCF fails, region=FALSE. If no variants, region=NULL
  print("This region does not contain any variants")

load(sprintf("gffseq_chr%s.RData",chrom)) # Annotation data from GFFtoFASTA

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
allwinsites <- slide@region.data@biallelic.sites # Positions in each window
allwindex <- slide@SLIDE.POS


#####################################
## GLOBAL STATISTICS AND VARIABLES ##
#####################################

## NEUTRALITY TESTS ##

# This calculates neutrality stats (Watterson's Theta, Fu and Li, Tajima's D) and the
# number of segregating sites (S). However, S excludes missing information in the outgroup.
slide <- neutrality.stats(slide)

## SITE FREQUENCY SPECTRA (IMPROVED) ##

n <- length(humans) # N = samples. Diploid, not divided by 2
bialhuman <- get.biallelic.matrix(region,1) [1:n,,drop=F] #  Whole biallelic matrix, remove outgroup
MAF <- colSums(bialhuman == 1)/nrow(bialhuman)
MAF.df <- data.frame(POS=names(MAF),MAF=unname(MAF)) # Position and MAF of biallelic variants

# IDENTIFY ANCESTRAL ALLELES (AA) FOR SNPS IN THE REGION
library(stringr)
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
mafan <- merge(MAF.df[,1,drop=F],mafan,all.x=TRUE) # Retain originals
# Frequencies do not exactly match those of the VCF file! Why? Rounding perhaps
# Some sites do not have AA. In some cases because they are divergent, but others are polymorphic:
# MAF.df[which(is.na(mafan$MAF))[!(which(is.na(mafan$MAF)) %in% which(MAF.df$MAF == 0)),]

######################################
## CUSTOM FUNCTION WINDOW BY WINDOW ##
######################################

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

  tabsum <- as.data.frame(matrix(numeric(nwin*15),ncol=15,nrow=nwin))
  colnames(tabsum) <- c("S","Pi","DAF","Divsites","D","K","Unknown","Alpha","Fisher","Psel","Pneu","Dsel","Dneu","msel","mneu")

  # LOOP TO READ THE WINDOWS IN SLIDE VARIABLE:
  for (window in 1:nrow(windows)) {
    # EXTRACT SEGREGATING SITES (EXCLUDING OUTGROUP)
    bial <- get.biallelic.matrix(slide,window) # Biallelic matrix of the window
    if (is.null(bial)||dim(bial)[2]==0) { # When no variants are detected
      newrow <- rep(0,15) # Empty rows
      tabsum[window,] <- newrow
      next
    }
    bialhuman <- bial[1:n,,drop=F] # Remove outgroup (drop = F to keep 1 dimension)
    misshuman <- colSums(is.na(bialhuman))>0 # Sites missing in humans (e.g. structural variants)
    poly.sites <- apply(bialhuman,2,sum)>0 & !misshuman # Sites polymorphic in humans w/o missing
    bialhuman <- bialhuman[,poly.sites,drop=FALSE] # Keep only polymorphic sites
    # DETERMINE M (EXCLUDING MISSING AND POLYALLELIC)
    winstart <- windows[window,1] # Select the start position in that window
    # Total number of sites: 1) Remove NA in humans 2) Remove polyallelic sites
    if (!length(region@region.data@polyallelic.sites) == 0) { # Make sure list exists
      multial <- region@region.data@polyallelic.sites[[1]] # Positions of all polyalleles
      multisites <- sum(!is.na(match(multial,winstart:(winstart+wsize)))) # N polyallelic sites in that window
    } else {multisites <- 0} # If not available, we assume 0
    m <- wsize-sum(misshuman,na.rm=T)-multisites
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

    if (is.numeric (mafan$DAF)) {
      winDAF <- as.vector(na.omit(mafan$DAF[windex]))
      DAF <- hist(winDAF,seq(0.0,1,0.05),plot=F)$counts
      DAF <- paste(DAF,collapse = ";")
    } else {DAF <- NA}

    ## DIVERGENCE ##

    # COUNT THE NUMBER OF UNKNOWNS IN THE OUTGROUP:
    unknowns <- sum(is.nan(bial[total,])) # Displayed as NaN in the biallelic matrix

    # COUNT THE NUMBER OF DIVERGENCE SITES:
    divtotal <- bial[total,] == 1 & !poly.sites # 1/1 in the outgroup and excluding polymorphic
    divsites <- sum (divtotal,na.rm=T) # Polymorphic cannot be counted for divergence
    # CALCULATE D AND K:
    mout <- m - unknowns
    D <- round(divsites/mout,7) # Observed divergence: Proportion of sites with divergent nucleotides
    K <- round(-3/4*log(1-4/3*D),7) # Real divergence: Jukes and Cantor model

    ## NATURAL SELECTION REGIMES ##

    winsites <- allwinsites[[window]] # Biallelic sites in this window
    bial.class <- gffseq[winsites] # GFF feature in each site

    ## MKT CALCULATION

    # STANDARD: Same result as MKT function if done with all syn.

    Pneu <- sum(poly.sites[bial.class==4],na.rm=TRUE) # Neutral=4-fold
    Psel <- sum(poly.sites[bial.class==0],na.rm=TRUE) # Selected=0-fold
    Dneu <- sum(divtotal[bial.class==4],na.rm=TRUE)
    Dsel <- sum(divtotal[bial.class==0],na.rm=TRUE)
    NI <- (Psel/Pneu)/(Dsel/Dneu)
    alpha <- 1-NI
    DoS <- Dsel/(Dsel+Dneu)-Psel/(Psel+Pneu) # Direction of selection
    contingency.std <- matrix(c(Psel,Pneu,Dsel,Dneu),c(2,2))
    test <- if(!is.na(sum(contingency.std))){
    if(sum(contingency.std)>0){fisher.test(contingency.std)$p.value}
    } else {NA}

    # ADAPTED:

    # 1) Use four-fold degenerated sites 2) Remove sites according to DAF

    # We use MAF for the adjustment because it is available for all positions, whereas some
    # do not have DAF due to the lack of AA. For DAF, please use winDAF.
    
    neuMAF <- as.vector(na.omit(MAF[bial.class==4 & poly.sites==TRUE]))
    selMAF <- as.vector(na.omit(MAF[bial.class==0 & poly.sites==TRUE]))
    
    Pneu.less5 <- sum(neuMAF<0.05) # P0 MAF less than 5%
    Pneu.more5 <- sum(neuMAF>0.05) # P0 MAF more than 5% 
    Psel.less5 <- sum(selMAF<0.05) # Pi MAF less than 5% 
    Psel.more5 <- sum(selMAF>0.05) # Pi MAF more than 5%
    
    Psel.neutral.less5 <- Psel*(Pneu.less5/Pneu) # Proportion of neutral within the MAF < 5% class
    Psel.neutral <- Psel.neutral.less5 + Psel.more5 # For alpha
    Psel.weak <- Psel.less5 - Psel.neutral.less5
    
    alpha.cor <- 1-(Psel.neutral/Pneu)*(Dneu/Dsel)
    contingency <- matrix(c(Psel.neutral,Pneu,Dsel,Dneu),c(2,2))
    test <- if(!is.na(sum(contingency))){
      if(sum(contingency)>0){fisher.test(contingency)$p.value}
    } else {NA}

    # OTHER ESTIMATORS:

    # To fully implement extended MKT, we need ms/mns, that is, the number of sites
    # of each class. Doing that would require modifying the 'set.synnonsyn2' function
    # to obtain all codons and check fold in every position (0,1,2)
    
    m.neu <- sum(bial.class == 4,na.rm=T) # neu=4-fold
    m.sel <- sum(bial.class == 0,na.rm=T) # sel=0-fold
    
    ## ADD NEW ROW ##
    newrow <- c(S,Pi(k,m,n),DAF,divsites,D,K,unknowns,alpha.cor,test,Psel.neutral,Pneu,Dsel,Dneu,m.sel,m.neu)
    print(newrow)
    tabsum[window,] <- newrow
  }
  return(tabsum)
}

## INTEGRATION OF NEUTRALITY, DIVERSITY AND DIVERGENCE METRICS ##

# Results from neutrality stats module are divided by number of sites
regiondata <- measures(slide)
S2 <- slide@n.segregating.sites[,1] # Segretaging sites excluding unknowns
Tajima_D <- round(slide@Tajima.D[,1]/wsize,7)
FuLi_F <- round(slide@Fu.Li.F[,1]/wsize,7)
theta <- round(slide@theta_Watterson[,1]/wsize,7)
if (exists("S2")) {regiondata <- cbind(regiondata[,1:2],theta,S2,Tajima_D,FuLi_F,regiondata[,3:ncol(regiondata)])
} else { regiondata <- cbind(regiondata[,1:2],theta=NA,S2=0,Tajima_D=0,FuLi_F=0,regiondata[,3:ncol(regiondata)]) }

######################
## DATA EXPORTATION ##
######################

# We add the chromosome info and the coordinates. For convenience in the next steps:
# a) Coordinates are converted to BED format: 0-based, subtract 1 from start column.
# b) Chromosome is expressed in "chrNN" format.
windows[,1] <- windows[,1]-1
export <- cbind(chr=rep(paste(c("chr",chrom),collapse=""),NROW(windows)),windows,regiondata)

write.table(export,file=sprintf("WindowData_chr%s.tab",chrom),quote=FALSE,sep="\t",row.names=F,append=TRUE,
col.names=!file.exists(sprintf("WindowData_chr%s.tab",chrom))) # Column names written if file does not exist

if("RMySQL" %in% rownames(installed.packages()) == TRUE) {

	suppressMessages(library(DBI))
	suppressMessages(library(RMySQL))
	
	con <- dbConnect(RMySQL::MySQL(),
        user="roger", password="RM333",
        dbname="PEGH", host="158.109.215.40")

	first <- !dbExistsTable(con,db)

	dbWriteTable(con,value=export,name=db,row.names=F,append=T)

	if (first == TRUE) # Remove if we want to concatenate various chromosomes
	dbSendQuery(con,sprintf("ALTER TABLE %s CHANGE COLUMN start start VARCHAR(30);",db))
	dbSendQuery(con,sprintf("ALTER TABLE %s ADD PRIMARY KEY (start);",db))

	on.exit(dbDisconnect(con))
}

print(Sys.time() - start.time)
