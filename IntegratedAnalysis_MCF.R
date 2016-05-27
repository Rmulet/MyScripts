#!/usr/bin/Rscript

# setwd("~/Documents/3_EpigenomicsData/Roadmap/Intraindividual")
# mode <- "Intraindividual"; group <- "STL003"
# setwd("~/Documents/3_EpigenomicsData/Roadmap/Interindividual")
# mode <- "Interindividual"; group <- "Adipose_Derived_Mesenchymal_Stem_Cell_Cultured_Cells"

# Combined Epigenetic Analysis Pipeline, v0.1 - This script combinesthe scripts'ChipSeqAnalysis.R' and 
# 'MethAnalysis.R' into a single one to perform integrated analysis of epigenome data. To do so, it 
# retrieves windows from the PEGH database ('Genomics' table) and calculates the variation in the 
# epigenomic landscape.Importantly, it assumes that the filenames follow the Roadmap pattern; files 
# from other sources will require modifications in the regex part.

# Nomenclature: 'Group' is the object of the study [donor/tissue]; 'sample' is each individual dataset in that group.

######################
## ARGUMENT PARSING ##
######################

args <- commandArgs(trailingOnly = TRUE) # Import arguments from command line

usage <- function() {
  cat ("\nProgram: Methylation Analysis Pipeline, v0.1\n")
  cat ("\nUsage: MethAnalysis.R [Mode] [Donor/Tissue]\n")
  cat ("File name format (Roadmap style): [University].[CellType].[Mark].[Donor]\n\n")
  quit()
}

if (length(args) < 1) {usage()}
if (args[1] == "-h" || args[1] == "--help") {usage()}

mode <- args[1] # Operation mode of the script
group <- args[2] # Donor/tissue group of the analysis

################################################
## ACCESS MYSQL DATABASE AND RETRIEVE WINDOWS ##
################################################

start <- Sys.time()
# WINDOWS IN MYSQL:
suppressMessages(library(RMySQL))
# Retrieve data from MySQL:
con <- dbConnect(RMySQL::MySQL(),
                 user="root", password="RM333",
                 dbname="PEGH", host="localhost")
res <- dbSendQuery(con, "SELECT chr,start,end FROM Genomics_Pilot")
windows <- dbFetch(res,n=-1) # n=-1 for no limit in records
invisible(dbClearResult(res)) # Frees resources associated with the query
windows[,2:3] <- apply(windows[,2:3],2,as.numeric) # Make the variable numeric
wsize <- windows[1,3]-windows[1,2] # Windows size
chr <- substr(windows[1,1],4,nchar(windows[1,1])) # Chromosome number
ntotal <- nrow(windows) # Number of windows

# BED is 0-based, but GRanges is 1-based.
if (suppressMessages(!require("GenomicRanges"))) {
  print ("The 'GenomicRanges' package is missing and will be installed")
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
  library("GenomicRanges")
}
gffwindows <- data.frame(chr=windows[,1],start=windows[,2]+1,end=windows[,3])
grwindows <- with(gffwindows,GRanges(chr,IRanges(start,end)))

#######################################
## FUNCTION TO CONVERT WIG TO BIGWIG ##
#######################################  

# WARNING: Requires the use of 'bwtool' and 'WigtoBigWig', downloadable here:
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
# https://github.com/CRG-Barcelona/bwtool/wiki
# For 'WigtoBigWig', chromosome sizes are needed. They can be retrieved with
# the aplication 'fetchChromSizes', available in the first directory.

wigtobw <- function() { # By default, takes all windows
    # IF NECESSARY, CREATE BIGWIG:
    filechr <- paste(substr(file,1,nchar(file)-2),chr,sep="")
    out <- paste(substr(file,1,nchar(file)-7),".bw",sep="")
      if (file.exists(out) == FALSE) {
        # Select only rows of the chromosome of interest [chr22]:
        system(sprintf("gunzip -c %s | awk '/chrom=chr%s/{p=1}/chrom=chr[^%s]/{p=0}p' > %s",file,chr,chr,filechr))
        # Transform wig to bigWig:
        system(sprintf("wigToBigWig %s hg19.chrom.sizes %s",filechr,out))
        file.remove(filechr)
        }
return(out)
}

######################################
## FUNCTION TO EXTRACT DATA FROM BW ##
######################################  

bwfraction <- function(ini=1,step=ntotal) {
  for (file in markfiles) {
    tissue.id <- str_match(file,pattern)[,2]
    donor.id <- str_match(file,pattern)[,5]
    if (anyNA(c(tissue.id,mark,donor.id))) {
      cat(sprintf("File %s does not follow the required pattern. Please check with '-h' or '--help'",file))
      next;
    }
    if (mode == "Interindividual" && group == tissue.id) {
      sample <- paste("ind_",donor.id,sep="")  
    } else if (mode == "Intraindividual" && group == donor.id) {
      sample <- abbreviate(str_replace_all(tissue.id,"_",""))
    } else {
      cat(sprintf("File %s not appropriate for the current analysis \n",file));next}
    ## EXTRACT DATA FROM BW ## 
    # Import as a table to R:
    system(sprintf("bwtool extract bed windows.bed %s %s -decimals=3",file,paste(file,".wn",sep="")))
    tab <- read.table(paste(file,".wn",sep=""),stringsAsFactors = FALSE)
    # WARNING: The table is imported to a variable called SAMPLE.
    suppressWarnings(assign(sample,lapply(strsplit(tab[,5],","),as.numeric),envir=.GlobalEnv))
    samples <- c(samples,sample)
    ## EXTRACT DATA FROM NARROWPEAK ##
    if (mark != "Bisulfite-Seq"){
      peakfile <- paste0(substr(file,1,nchar(file)-16),"narrowPeak.gz")
      temp <- system(sprintf("gunzip -c %s | grep '%s' | cut -f1,2,3 |  sort",peakfile,windows[1,1]),intern=TRUE) 
      peaks <- read.table(textConnection(temp),sep="\t")
      names(peaks) <- c("chr","start","end")
      assign(paste0(sample,".peak"),with(peaks,GRanges(chr,IRanges(start,end))),envir=.GlobalEnv) # Transform to IRanges format
    }
  }
  return(samples)
}
# WE CAN INTEGRATE THIS IF THERE IS ONLY ONE SEGMENT OF WINDOWS   

################################################
## FUNCTION TO MEASURE EPIGENETIC VARIABILITY ##
################################################

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

episcore <- function(nwin=ntotal) { # By default, it takes all the windows
  # Preallocate vectors containing the results
  meanavg <- numeric(nwin); varavg <- numeric(nwin); meanvar <- numeric(nwin)
  pi <- numeric(nwin); theta <- numeric(nwin); nsites <- numeric(nwin)
  # Compares the n-line of each window across samples (tissues/donors/species)
  for (window in 1:nwin) {
    mat <- matrix(nrow=0,ncol=wsize) # Epigenetic diversity matrix
    matpeaks <- matrix(nrow=0,ncol=wsize)  # For Chip-Seq only
    for (sample in samples) { # We add every sample to the matrix
      # For Chip-Seq and Methylation, we retrieve intensity data
      mat <- rbind(mat,eval(as.symbol(sample))[[window]])
      # For Chip-Seq, we identify with 1 positions of the window that have a peak
      if (mark != "Bisulfite-Seq") {
        newrow <- numeric(wsize)
        grpeaks <- eval(as.symbol(paste0(sample,".peak")))
        overlap <- findOverlaps(query=grwindows[window],subject=grpeaks)
        spans <- ranges(overlap,ranges(grwindows[window]),ranges(grpeaks)) # Ranges of the overlap
        if (length(spans) > 0) {
          for (i in 1:length(spans)) { # For each peak overlapping that window
            newrow[(start(spans[i])-start(grwindows[window])):(end(spans[i])-start(grwindows[window]))] <- 1
          }
        }  
        matpeaks <- rbind(matpeaks,newrow)
        rownames(matpeaks)[nrow(matpeaks)] <- sample
      }
    }
    matclean <- mat[,complete.cases(t(mat))] # We remove positions with NA (non-meth/missing)
    ## VARIANCE AND SIGNAL INTENSITY ## 
    if (dim(matclean)[2] == 0) { # No methylated positions in this region
      meanavg[window] <- 0
      varavg[window] <- 0
      meanvar[window] <- 0
      next
    }
    avg <- apply(matclean,1,mean) # Average of methylation in METHYLATION LOCI
    meanavg[window] <- mean(avg) # Average of all samples
    varavg[window] <- var(avg) # Variance of average across samples
    posvar <- apply(matclean,2,var) # Variance across samples in each position
    meanvar[window] <- mean(posvar) # Average variance across samples
    ## METHYLATION POLYMORPHISM ESTIMATORS ##
    if (mark == "Bisulfite-Seq") {
      episites <- apply(matclean,c(1,2),function(x){if(x>=0.7){x <- 1} # Assigns 1/0/0.5
        else if(x < 0.7 && x > 0.3){x <- 0.5} else {0} })
      m <- ncol(cpgsites) # Total CpG sites in the window (methylation loci)
      n <- length(samples)*2 # "Diploid"
      # ECF = Epigenetic Call Format: Only variable positions
      ECF <- cpgsites[,apply(cpgsites,2,function(x){!all(x==x[1])}),drop=FALSE]
      ki <- apply(ECF,2,function(x) {a <- table(x);
       (sum(a[names(a) == 0])*2+sum(a[names(a) == 0.5]))*
       (sum(a[names(a) == 1])*2+sum(a[names(a) == 0.5]))})
      ## OTHER EPIGENETIC MARKS POLYMORPHISM ESTIMATORS ##
    } else {
      m <- ncol(matpeaks) # Total nucleotides in that window (1000)
      n <- nrow(matpeaks) # "Haploid" n
      ECF <- matpeaks[,apply(matpeaks,2,function(x){!all(x==x[1])}),drop=FALSE] # Can be joined with MCF
      ki <- apply(ECF,2,table)
      print(ECF)
    }
    if (is.null(bialhuman)||dim(ECF)[2] == 0) { # No differentially methylated positions in the region
      S <- 0
      ki <- 0
    } else {  
      k <- sum(ki)
      nsites <- ncol(ECF) # Segregating sites = positions with epigenetic differences
      print(nsites)
      pi[window] <- Pi(k,m,n)
      theta[window] <- Theta(S,m,n)
      nsites[window] <- nsites
    }
 }
 epidata <- data.frame(Pi=pi,Theta=theta,S=nsites,Level=meanavg,LevelVar=varavg,Var=meanvar)
 return(epidata)
}

#############################
## ANALYZE FOLDER CONTENTS ##
#############################

library(stringr)

filenames <- list.files(".", pattern=sprintf("%s.+(bigwig$|wig\\.gz$)",group), full.names=TRUE) # Files in the folder
# Roadmap standard: a few donors have IDs with dots 
pattern <- "\\.(\\w+)\\.(Bisulfite-Seq|H[A|2B|3|4]K\\d+(me\\d|ac))\\.(\\w+)\\."
samples <- list()
marks <- vector()

# GROUP THE FILES BY EPIGENETIC MARK:
for (file in filenames) {
  mark <- str_match(file,pattern)[,3]
  if (is.na(mark)==TRUE) {
    cat(sprintf("File %s does not follow the required pattern. Please check with '-h' or '--help' \n",file))
    next;
  }
  if (!mark %in% marks) {
    marks <- c(marks,mark)
  }
  if (mark == "Bisulfite-Seq") {
    wigtobw()
  }
  filenames <- list.files(".", pattern=sprintf("%s.+(bigwig$|bw$)",group), full.names=TRUE) # BigWig files
}

# CALCULATE VARIABILITY FOR EACH MARK:
epidata <- windows
for (mark in marks){
  markfiles <- grep(mark,filenames,value=T)
  samples <- vector()
  # Verifying the number of samples for mark (can be done above):
  if (length(markfiles) < 2) {
    if (mode == "Intraindividual") {
      stop (sprintf("Please provide %s data from more than one tissue",mark))
    } else if (mode == "Interindividual") {
      stop (sprintf("Please provide %s data from more than one donor",mark))
    }
    # FRACTIONED ANALYSIS:  
  } else if (ntotal*wsize < 500000) { # If the region contains < 500 kb, all is processed at once
    write.table(windows,file="windows.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    samples <- bwfraction()
    markdata <- episcore()
  } else { # If the region contains > 500 kb, the analysis is fractioned
    markdata <- data.frame(Level=numeric(0),LevelVar=numeric(0),Var=numeric(0))
    chunk <- 500000/wsize
    for (ini in seq(1,ntotal,by=chunk)) { # REVIEW INTERVAL!!!!
      if ((n-ini+1) < chunk) {
        # Create 'windows.bed' file for 'bwtools':
        chunk <- ntotal%%chunk
      }
      write.table(windows[ini:(ini+chunk-1),],file="windows.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
      bwfraction(ini=ini,step=chunk)
      markdata <- rbind(markdata,episcore(nwin=chunk))
    }
  }
  labels <- c(sprintf("%s_Level",mark),sprintf("%s_LevelVar",mark),sprintf("%s_Var",mark))
  colnames(markdata) <- labels
  epidata <- cbind(epidata,markdata) # If no fractions, that can be moved above
}

if (mode == "Intraindividual") {
  code <- group  # Samples of the study
} else if (mode == "Interindividual") {
  code <- abbreviate(str_replace_all(tissue,"_",""))
} else if (mode == "Interspecies") {
  code <- abbreviate(str_replace_all(tissue,"_",""))}

#############################
## DATA EXPORT TO MYSQL DB ##
#############################

marks <- paste(marks,collapse=";")
epidata[,4:ncol(epidata)] <- apply(epidata[,4:ncol(epidata)],2,round,5) # ROUND(OPTIONAL)
epidata.name <- paste(substr(mode,1,7),"_",code,sep="")

dbWriteTable(con,epidata.name,epidata,row.names=F,overwrite=T)

## SUMMARY TABLE ##

if (mode == "Interindividual") {
  donorlist <- unique(samples) # Combine samples from chip-seq and methylation
  sqlquery <- paste("INSERT INTO Interindividual VALUES('"
                    ,epidata.name,"','",code,"','",NROW(samples),"','hg19','Roadmap','",marks,"');",sep="")
  print(samples)
} else if (mode == "Intraindividual") {
  # Combine samples from chip-seq and methylation
  samples <- data.frame(Tissue=names(samples),Abbreviation=unname(samples))
  samples <- unique(samples)
  sqlquery <- paste("INSERT INTO Intraindividual VALUES('"
                    ,epidata.name,"','",code,"','","Unknown","','",NROW(samples),"','hg19','Roadmap','",marks,"');",sep="")
  print(samples)
  dbWriteTable(con,"TissueAbbreviations",samples,append=T,row.names=F)
  dbSendQuery(con,"ALTER TABLE TissueAbbreviations ADD PRIMARY KEY (Tissue);")
}
invisible(dbSendQuery(con,sqlquery))
invisible(dbDisconnect(con))

Sys.time() - start
cat(c(ntotal*wsize,"bases"))


#  The -log10(p-value) scores provide a convenient way to threshold signal (e.g. 2 corresponds to a p-value threshold of 1e-2), similar to what is used in identifying enriched regions (peak calling). We recommend using the signal confidence score tracks for visualization. A universal threshold of 2 provides good separation between signal and noise. Both types of signal tracks were also generated for the unconsolidated datasets using the same parameter settings described above.
 