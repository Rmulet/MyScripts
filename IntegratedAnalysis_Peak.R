#!/usr/bin/Rscript

# setwd("~/Documents/3_EpigenomicsData/Roadmap/Intraindividual/NarrowPeak")
# mode <- "Intraindividual"; group <- "STL003"
# setwd("~/Documents/3_EpigenomicsData/Roadmap/Interindividual/NarrowPeak")
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
res <- dbSendQuery(con, "SELECT chr,start,end FROM Genomics")
windows <- dbFetch(res)
invisible(dbClearResult(res)) # Frees resources associated with the query
windows[,2:3] <- apply(windows[,2:3],2,as.numeric) # Make the variable numeric
wsize <- windows[1,3]-windows[1,2] # Windows size
chr <- substr(windows[1,1],4,nchar(windows[1,1])) # Chromosome number
ntotal <- nrow(windows) # Number of windows

# GRANGES group (FOR CHIP-SEQ)
if (suppressMessages(!require("GenomicRanges"))) {
  print ("The 'GenomicRanges' package is missing and will be installed")
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
  library("GenomicRanges")
}
# BED is 0-based, but GRanges is 1-based
gffwindows <- data.frame(chr=windows[,1],start=windows[,2]+1,end=windows[,3])
grwindows <- with(gffwindows,GRanges(chr,IRanges(start,end)))

#####################################################
## FIND OVERLAPS BETWEEN WINDOWS AND CHIP-SEQ DATA ##
#####################################################

chipseqscore <- function(file) {
  options(digits=3)
  # IMPORT AND PROCESS THE FILE
  # Processing the file with system commands is much faster than R.
  temp <- system(sprintf("gunzip -c %s | grep '%s' | cut -f1,2,3 |  sort",file,windows[1,1]),intern=TRUE) 
  peaks <- read.table(textConnection(temp),sep="\t")
  names(peaks) <- c("chr","start","end")
  grpeaks <- with(peaks,GRanges(chr,IRanges(start,end))) # Transform to IRanges format
  
  # IDENTIFY OVERLAPS BETWEEN PEAKS AND WINDOWS
  overlap <- findOverlaps(query=grwindows,subject=grpeaks)
  spans <- width(ranges(overlap,ranges(grwindows),ranges(grpeaks))) # Measures span of the overlap
  query <- queryHits(overlap)
  
  # SOLVE PEAK OVERLAP IN THE SAME WINDOW
  qdup <- unique(query[duplicated(query)]) # Windows with more than one peak
  for (i in qdup) { 
    score <- 0
    rep <- which(query == i)
    spanr <- spans[rep]
    peakspos <- subjectHits(overlap)[rep]
    tab <- peaks[peakspos,]
    for (n in 2:length(rep)) {
      if (tab[n,2] > max(tab[1:n-1,3])) { # Non-overlapping peaks
        score <- spanr[n-1]+spanr[n]
      }
      if (tab[n,2] < max(tab[1:n-1,3])) { # Overlapping peaks
        if (tab[n,3] < max(tab[1:n-1,3])) { # New peak contained the former
          next
        }
        else if (tab[n,3] > max(tab[1:n-1,3])) { # New peak extends beyond the former
          score <- spanr[n-1]+(spanr[n]-(tab[n,3]-max(tab[1:n-1,3])))
        }
      }
    }
    # 
    query <- query[-rep[2:length(rep)]] # Remove all repeated occurrences of the window
    spans[rep[1]] <- score # Replace the span value of the first occurrence of the repeat
    spans <- spans[-rep[2:length(rep)]] # Remove the values of the repeated occurrences
  }  
  # PRINT VALUES 
  values <- spans/wsize*100 # Percentage of overlap in the interval
  percentages <- numeric(queryLength(overlap))
  percentages[query] <- values
  return(percentages)
}

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
        }
return(out)
}

# EXTRACT INFORMATION IN WINDOWS:
bwfraction <- function(step=ntotal) {
  for (i in 1:length(metfiles)) {
    file <- metfiles[i]; sample <- metsamples[i]
    # Import as a table to R:
    system(sprintf("bwtool extract bed windows.bed %s %s -decimals=3",file,paste(file,".wn",sep="")))
    tab <- read.table(paste(file,".wn",sep=""),stringsAsFactors = FALSE)
    # WARNING: The table is imported to a variable called SAMPLE.
    suppressWarnings(assign(sample,lapply(strsplit(tab[,5],","),as.numeric),envir=.GlobalEnv))
  }
}
# WE CAN INTEGRATE THIS IF THERE IS ONLY ONE SEGMENT OF WINDOWS   

############################################
## FUNCTION TO MEASURE WINDOW METHYLATION ##
############################################

methylscore <- function(nwin=ntotal) { # By default, it takes all the windows
  # Preallocate vectors containing the results
  meanavg <- numeric(nwin)
  varavg <- numeric(nwin)
  meanvar <- numeric(nwin)
  # Compares the n-line of each window across samples (tissues/donors/species)
  for (window in 1:nwin) {
    mat <- matrix(nrow=0,ncol=wsize) # Epigenetic diversity matrix
    for (sample in metsamples) { # We add every sample to the matrix
      mat <- rbind(mat,eval(as.symbol(sample))[[window]])
    }
    matmet <- mat[,complete.cases(t(mat))] # We remove NA (non-meth/missing)
    if (dim(matmet)[2] == 0) { # No methylated positions in this region
      meanavg[window] <- 0
      varavg[window] <- 0
      meanvar[window] <- 0
      next
    }
    avg <- apply(matmet,1,mean) # Average of methylation in METHYLATION LOCI
    meanavg[window] <- mean(avg)
    varavg[window] <- var(avg)
    posvar <- apply(matmet,2,var)
    meanvar[window] <- mean(posvar)
  }
  methyldata <- data.frame(MethLevel=meanavg,MethLevelVar=varavg,MethVar=meanvar)
  return(methyldata)
}

#############################
## ANALYZE FOLDER CONTENTS ##
#############################

library(stringr)

filenames <- list.files(".", pattern=sprintf("%s.+(narrowPeak|wig\\.gz)",group), full.names=TRUE) # Files in the folder
# Roadmap standard: a few donors have IDs with dots 
pattern <- "\\.(\\w+)\\.(Bisulfite-Seq|H[A|2B|3|4]K\\d+(me\\d|ac))\\.(\\w+)\\."

if (mode == "Intraindividual") {
  code <- group  # Samples of the study
} else if (mode == "Interindividual") {
  code <- abbreviate(str_replace_all(tissue,"_",""))
} else if (mode == "Interspecies") {
  code <- abbreviate(str_replace_all(tissue,"_",""))}

tablist <- vector()
samples <- vector()
metfiles <- vector()
metsamples <- vector()

for (file in filenames) {
  tissue.id <- str_match(file,pattern)[,2]
  mark <- str_match(file,pattern)[,3]
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
  if (mark == "Bisulfite-Seq") {
    print(sample)
    metfiles <- c(metfiles,wigtobw())
    metsamples <- c(metsamples,sample)
    } else {
    samples <- c(samples,sample)
    tab.id <- paste(c("intrain_",code,"_",mark),collapse="") # It contains the NAME of the variable
    if ((tab.id %in% tablist)==FALSE) {
      tablist <- c(tablist,tab.id)       
      assign(tab.id,windows)
    }
    tab.res <- eval(as.symbol(tab.id))
    tab.res <- data.frame(tab.res,chipseqscore(file))
    colnames(tab.res)[ncol(tab.res)] <- sample
    assign(tab.id,tab.res)
  }
}

###############################
## CHIP-SEQ HISTONE ANALYSIS ##
###############################

chipdata <- windows
for (tab in tablist) { 
  tab.res <- eval(as.symbol(tab))
  if (ncol(tab.res) < 5) { # Histones analyzed in only one tissue/individual (can be adjusted)
    if (mode == "Intraindividual") {
      cat (sprintf("Please provide more than one tissue for histone %s\n",unlist(str_split(tab,"\\_"))[3]))
    } else if (mode == "Interindividual") {
      cat (sprintf("Please provide more than one donor for histone %s\n",unlist(str_split(tab,"\\_"))[3]))
    }
    next
  }
  # We upload one table for each histone and calculate the mean and the variance
  # dbWriteTable(con,name=tab,value=eval(as.symbol(tab)),row.names=F,overwrite=T)
  markmean <- apply(tab.res[,-(1:3)],1,mean)
  markvar <-  apply(tab.res[,-(1:3)],1,var)
  chipdata <- cbind(chipdata,markmean,markvar)
  colnames(chipdata)[ncol(chipdata)-1] <- paste(unlist(str_split(tab,"\\_"))[3],"_mean",sep="")
  colnames(chipdata)[ncol(chipdata)] <- paste(unlist(str_split(tab,"\\_"))[3],"_var",sep="")
}
# IF NO CHIP-SEQ FILES, THEN CHIPDATA = WINDOWS

##########################
## METHYLATION ANALYSIS ##
##########################

# Verifying the number of samples for mark (can be done above):
if (length(metfiles) < 2) {
  if (mode == "Intraindividual") {
    stop ("Please provide methylation data from more than one tissue")
  } else if (mode == "Interindividual") {
    stop ("Please provide methylation data from more than one donor")
  }
# Determining whether the analysis should be fractioned:  
} else if (ntotal*wsize < 500000) { # If the region contains < 500 kb, all is processed at once
  write.table(windows,file="windows.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  bwfraction()
  methyldata <- methylscore()
} else { # If the region contains > 500 kb, the analysis is fractioned
  methyldata <- data.frame(MethLevel=numeric(0),MethLevelVar=numeric(0),MethVar=numeric(0))
  chunk <- 100
  for (ini in seq(1,ntotal,by=chunk)) { # REVIEW INTERVAL!!!!
    if ((n-ini+1) < chunk) {
      # Create 'windows.bed' file for 'bwtools':
      chunk <- ntotal%%chunk
    }
    write.table(windows[ini:(ini+chunk-1),],file="windows.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    bwfraction(step=chunk)
    methyldata <- rbind(methyldata,methylscore(nwin=chunk))
  }
}
Sys.time() - start

#############################
## DATA EXPORT TO MYSQL DB ##
#############################

## RESULTS TABLE (TISSUE/DONOR) ##
marks <- unique(sapply(strsplit(colnames(chipdata)[4:length(colnames(chipdata))],"_"),'[',1))
# If no methylation or chipseq, the 'result' table should also contain the existing data
if (exists("methyldata") == TRUE) {
  marks <- c(marks,"Methylation")
  result <- data.frame(chipdata,methyldata)
} else {result <- chipdata}
marks <- paste(marks,collapse=";")
result[,4:ncol(result)] <- apply(result[,4:ncol(result)],2,round,5) # ROUND(OPTIONAL)
result.name <- paste(substr(mode,1,7),"_",code,sep="")

dbWriteTable(con,result.name,result,row.names=F,overwrite=T)

## SUMMARY TABLE ##

if (mode == "Interindividual") {
  donorlist <- unique(c(samples,metsamples)) # Combine samples from chip-seq and methylation
  sqlquery <- paste("INSERT INTO Interindividual VALUES('"
                    ,result.name,"','",code,"','",NROW(samples),"','hg19','Roadmap','",marks,"');",sep="")
  print(samples)
} else if (mode == "Intraindividual") {
  # Combine samples from chip-seq and methylation
  samples <- data.frame(Tissue=c(names(samples),names(metsamples)),Abbreviation=c(unname(samples),unname(metsamples)))
  samples <- unique(samples)
  sqlquery <- paste("INSERT INTO Intraindividual VALUES('"
                    ,result.name,"','",code,"','","Unknown","','",NROW(samples),"','hg19','Roadmap','",marks,"');",sep="")
  print(samples)
  dbWriteTable(con,"TissueAbbreviations",samples,append=T,row.names=F)
}
invisible(dbSendQuery(con,sqlquery))

invisible(dbDisconnect(con))
Sys.time() - start
cat(c(ntotal*wsize,"bases"))