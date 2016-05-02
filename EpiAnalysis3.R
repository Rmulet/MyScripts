#!/usr/bin/Rscript

# setwd("~/Documents/3_EpigenomicsData/Roadmap")
# file <- "BI.Adipose_Nuclei.H3K4me1.7.narrowPeak"

# Epigenetic Analysis Pipeline, v0.3 - This script heavly relies on the use of BED files
# for narrow peaks, so files must be provided in the specified format. Importantly,
# iit assumes that the filenames follow the Roadmap pattern; files from other sources
# will require modifications in the regex part.

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
  options(digits=3)
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
# Roadmap standard: a few donors have IDs with dots 
pattern <- "\\.(\\w+)\\.(H[A|2B|3|4]K\\d(me\\d|ac))\\.(.+)\\.narrowPeak" 
tablist <- vector()

for (file in filenames) {
  tissue.id <- str_match(file,pattern)[,2]
  mark <- str_match(file,pattern)[,3]
  donor.id <- str_match(file,pattern)[,5]
  if (mode == "Interindividual" && tissue == tissue.id) {
    code <- abbreviate(str_replace_all(tissue,"_",""))
    tab.id <- paste(c("interin.",code,".",donor.id),collapse="") # It contains the NAME of the variable
    if ((tab.id %in% tablist)==FALSE) {
      tablist <- c(tablist,tab.id)       
      assign(tab.id,windows)
    }
    tab.res <- eval(as.symbol(tab.id))
    tab.res <- data.frame(tab.res,chipseqscore(file))
    colnames(tab.res)[ncol(tab.res)] <- mark
    assign(tab.id,tab.res)
  }
}

###################################
## EPIGENOMIC DIVERSITY ANALYSIS ##
###################################

for (sample in tablist) {
  
}


######################
## FEATURE ANALYSIS ##
######################

# This module of the function attempts to find a correlation between windows with
# epigenetic marks and genomic features such as genes. 

system("cat chr22.gff | sed 's/^22/chr22/g' > Chr22.gff")
gff <- read.table('Chr22.gff',header=FALSE)[,c(1,3,4,5)]
names(gff) <- c("chr","feature","start","end")

evaluate <- function(gff,selection="gene",mode="any") { 

grgff <- with(gff,GRanges(chr,IRanges(start,end),feature=feature))
grselect <- grgff[grgff$feature==selection]
overlap <- findOverlaps(query=grwindows,subject=grselect)
  
feature <- rep(0,queryLength(overlap))
feature[unique(queryHits(overlap))] <- selection
richmarks <- cbind(result,feature)
  
if (mode=="any") { # Considers marks on ANY histone
  
  # Windows with peaks in at least one mark + SELECTION
  subset(richmarks,(H3K4me1>0|H3K4me1.1>0|H3K4me1.2>0|H3K9ac>0)&feature==selection)
  a <- nrow(subset(richmarks,(H3K4me1>0|H3K4me1.1>0|H3K4me1.2>0|H3K9ac>0)&feature==selection))
  # Windows with peaks in at least one mark - SELECTION
  b <- nrow(subset(richmarks,(H3K4me1>0|H3K4me1.1>0|H3K4me1.2>0|H3K9ac>0)&feature=="0"))
  # Windows without marks + SELECTION:
  c <- nrow(subset(richmarks,(H3K4me1==0&H3K4me1.1==0&H3K4me1.2==0&H3K9ac==0)&feature==selection))
  # Windows without marks - SELECTION:
  d <- nrow(subset(richmarks,(H3K4me1==0&H3K4me1.1==0&H3K4me1.2==0&H3K9ac==0)&feature=="0"))
  }

if (mode=="all") {
  # Windows with peaks in all marks + SELECTION
  a <- nrow(subset(richmarks,H3K4me1>0&H3K4me1.1>0&H3K4me1.2>0&H3K9ac>0&feature==selection))
  # Windows with peaks in all marks - SELECTION
  b <- nrow(subset(richmarks,H3K4me1>0&H3K4me1.1>0&H3K4me1.2>0&H3K9ac>0&feature=="0"))
  # Windows without marks + SELECTION:
  c <- nrow(subset(richmarks,(H3K4me1==0&H3K4me1.1==0&H3K4me1.2==0&H3K9ac==0)&feature==selection))
  # Windows without marks - SELECTION:
  d <- nrow(subset(richmarks,(H3K4me1==0&H3K4me1.1==0&H3K4me1.2==0&H3K9ac==0)&feature=="0"))
  }
  
# STATISTICAL TEST:
  evaluation <- matrix(nrow=2,ncol=2,c(a,b,c,d),dimnames=list(c("Feat","NoFeat"),c("Mark","NoMark")))
  print(evaluation)
  if (any(c(a,b,c,d)<5)) {
    print(fisher.test(evaluation)) # P-value is 1
  } else {
    print(chisq.test(evaluation))
  }
}
  
invisible(dbDisconnect(con))