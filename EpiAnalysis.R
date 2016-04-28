#!/usr/bin/Rscript

# setwd("~/Documents/3_EpigenomicsData/Roadmap")
# file <- "BI.Adipose_Nuclei.H3K4me1.7.narrowPeak"

# Epigenetic Analysis Pipeline, v0.3 - This script heavly relies on the use of BED files
# for narrow peaks, so files must be provided in the specified format.

# WARNING: Review positions!

######################
## ARGUMENT PARSING ##
######################

args <- commandArgs() # Import arguments from command line
mode <- args[6] # Operation mode of the script
if (mode == "Intraindividual") {
  donor <- args[7] # Name of the donor
  print("Intraindividual mode selected")
} else if (mode == "Interindividual") {
  tissue <- args[7] 
  print("Interindividual mode selected")
} else if (mode == "Interspecies") {
  tissue <- args[7]
  print("Interspecies mode selected")
}

headerbed <- c("chr","start","end")

################################################
## ACCESS MYSQL DATABASE AND RETRIEVE WINDOWS ##
################################################

if (suppressMessages(!require("GenomicRanges"))) {
  print ("The 'GenomicRanges' package is missing and will be installed")
  source("https://bioconductor.org/biocLite.R")
  biocLite("GenomicRanges")
}

# WINDOWS IN MYSQL:
suppressMessages(library(RMySQL))
# Retrieve data from MySQL and put in GRanges format
con <- dbConnect(RMySQL::MySQL(),
                 user="root", password="RM333",
                 dbname="PEGH", host="localhost")
res <- dbSendQuery(con, "SELECT Chr,Window FROM Genomics")
windows <- dbFetch(res)
invisible(dbClearResult(res)) # Frees resources associated with the query
coords <- trimws(do.call(rbind,strsplit(windows[,2],"-"))) # Separate the coordinates bound by "-"
coords <- apply(coords,2,as.numeric) # Make the variable numeric
windows <- data.frame(paste("chr",windows[,1],sep=""),coords)
names(windows) <- headerbed
grwindows <- with(windows,GRanges(chr,IRanges(start,end)))

#####################################################
## FIND OVERLAPS BETWEEN WINDOWS AND CHIP-SEQ DATA ##
#####################################################

chipseqscore <- function(file) {
  
  # IMPORT AND PROCESS THE FILE
  # Processing the file with system commands is much faster than R.
  temp <- system(sprintf("gunzip -c %s | grep 'chr22' | cut -f1,2,3 |  sort",file),intern=TRUE) 
  peaks <- read.table(textConnection(temp),sep="\t")
  names(peaks) <- headerbed
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
  values <- spans/200*100
  percentages <- numeric(queryLength(overlap))
  percentages[query] <- values
  return(percentages)
}

#############################
## ANALYZE FOLDER CONTENTS ##
#############################

library(stringr)

filenames <- list.files(".", pattern="*.narrowPeak", full.names=TRUE) # Files in the folder
pattern <- "\\.(\\w+)\\.(H[A|2B|3|4]K\\d(me\\d|ac))\\.(\\w+)" # Roadmap standard
result <- windows

for (file in filenames) {
  tissue.id <- str_match(file,pattern)[,2]
  mark <- str_match(file,pattern)[,3]
  donor.id <- str_match(file,pattern)[,5]
  
  if (mode == "Interindividual" && tissue == tissue.id) {
    assign(paste(c("tab.",donor.id),collapse=""),windows)
    result <- data.frame(result,chipseqscore(file))
    colnames(result)[ncol(result)] <- mark
  }
}

print(head(result))
invisible(dbDisconnect(con))