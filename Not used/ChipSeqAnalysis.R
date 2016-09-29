#!/usr/bin/Rscript

# setwd("~/Documents/3_EpigenomicsData/Roadmap/Intraindividual")
# mode <- "Intraindividual"; donor <- "STL003"

# Epigenetic Analysis Pipeline, v0.3 - This script compares windows contained in a MySQL database with
# ChIP-seq data in narrowPeak format to investigate what fraction of each window falls within a peak.
# Importantly, it assumes that the filenames follow the Roadmap pattern; files from other sources
# will require modifications in the regex part.

# WARNING: Review positions!

######################
## ARGUMENT PARSING ##
######################

args <- commandArgs(trailingOnly = TRUE) # Import arguments from command line

usage <- function() {
  cat ("\nProgram: Epigenetic Analysis Pipeline, v0.3 - This script compares windows contained in a MySQL database with ChIP-seq data in narrowPeak format to investigate what fraction of each window falls within a peak.\n")
  cat ("\nUsage: EpiAnalysis.R [Mode] [Donor/Tissue]\n")
  cat ("File name format (Roadmap style): [University].[CellType].Bisulfite-Seq.[Donor]\n\n")
  quit()
}

if (length(args) < 1) {print("Please provide an adequate number of arguments");usage()}
if (args[1] == "-h" || args[1] == "--help") {usage()}

mode <- args[1] # Operation mode of the script
if (mode == "Intraindividual") {
  donor <- args[2] # Name of the donor
  cat("Intraindividual mode selected")
} else if (mode == "Interindividual") {
  tissue <- args[2] 
  cat("Interindividual mode selected")
} else if (mode == "Interspecies") {
  tissue <- args[2]
  cat("Interspecies mode selected")
}

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
res <- dbSendQuery(con, "SELECT chr,start,end FROM Genomics")
windows <- dbFetch(res)
invisible(dbClearResult(res)) # Frees resources associated with the query
windows[,2:3] <- apply(windows[,2:3],2,as.numeric) # Make the variable numeric
wsize <- windows[1,3]-windows[1,2] # Windows length

grwindows <- with(windows,GRanges(chr,IRanges(start,end)))

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

#############################
## ANALYZE FOLDER CONTENTS ##
#############################

library(stringr)

filenames <- list.files(".", pattern="*.narrowPeak", full.names=TRUE) # Files in the folder
# Roadmap standard: a few donors have IDs with dots 
pattern <- "\\.(\\w+)\\.(H[A|2B|3|4]K\\d+(me\\d|ac))\\.(\\w+)\\.narrowPeak" 
tablist <- vector()
donorlist <- vector()
tissuelist <- data.frame(Tissue=character(0),Abbreviation=character(0),stringsAsFactors=FALSE)

for (file in filenames) {
  tissue.id <- str_match(file,pattern)[,2]
  mark <- str_match(file,pattern)[,3]
  donor.id <- str_match(file,pattern)[,5]
  if (anyNA(c(tissue.id,mark,donor.id))) {
    stop(sprintf("File %s does not follow the required pattern. Please check with '-h' or '--help'",file))
  }
  if (mode == "Interindividual" && tissue == tissue.id) {
    code <- abbreviate(str_replace_all(tissue.id,"_","")) # Summary table
    tab.id <- paste(c("interin_",code,"_",mark),collapse="") # It contains the NAME of the variable
    col.id <- paste("ind_",donor.id,sep="")  
    donorlist <- c(donorlist,donor.id)
  } else if (mode == "Intraindividual" && donor == donor.id) {
    code <- donor.id  # Summary table
    tab.id <- paste(c("intrain_",code,"_",mark),collapse="") # It contains the NAME of the variable
    col.id <- abbreviate(str_replace_all(tissue.id,"_",""))
    tissuelist[nrow(tissuelist)+1,] <- c(tissue.id,col.id)
  } else if (mode == "Interspecies" && tissue == tissue.id) {
    code <- abbreviate(str_replace_all(tissue.id,"_","")) # Summary table
    tab.id <- paste(c("intersp_",code,"_",mark),collapse="") # It contains the NAME of the variable
    col.id <- paste("ind_",donor.id,sep="")  
    donorlist <- c(donorlist,donor.id)
  } else {next} # If it does not fulfill any condition above, next cycle
  
  if ((tab.id %in% tablist)==FALSE) {
    tablist <- c(tablist,tab.id)       
    assign(tab.id,windows)
  }
  tab.res <- eval(as.symbol(tab.id))
  tab.res <- data.frame(tab.res,chipseqscore(file))
  colnames(tab.res)[ncol(tab.res)] <- col.id
  assign(tab.id,tab.res)
}

print(tablist)
donorlist <- unique(donorlist)
tissuelist <- unique(tissuelist)

###################################
## EPIGENOMIC DIVERSITY ANALYSIS ##
###################################

result <- windows
for (tab in tablist) { 
  tab.res <- eval(as.symbol(tab))
  if (length(tab.res) < 5) { # Histones analyzed in only one tissue/individual (can be adjusted)
    if (mode == "Intraindividual") {
      cat (sprintf("Please provide more than one tissue for histone %s\n",unlist(str_split(tab,"\\_"))[3]))
    } else if (mode == "Interindividual") {
      cat (sprintf("Please provide more than one donor for histone %s\n",unlist(str_split(tab,"\\_"))[3]))
    }
   next
  }
  # We upload one table for each histone and calculate the mean and the variance
  dbWriteTable(con,name=tab,value=eval(as.symbol(tab)),row.names=F,overwrite=T)
  markmean <- apply(tab.res[,-(1:3)],1,mean)
  markvar <-  apply(tab.res[,-(1:3)],1,var)
  result <- cbind(result,markmean,markvar)
  colnames(result)[ncol(result)-1] <- paste(unlist(str_split(tab,"\\_"))[3],"_mean",sep="")
  colnames(result)[ncol(result)] <- paste(unlist(str_split(tab,"\\_"))[3],"_var",sep="")
}
# Then we upload the mean and the variance in that histone
result[,4:ncol(result)] <- apply(result[,4:ncol(result)],2,round,5) # WARNING: OPTIONAL
result.name <- paste(substr(mode,1,7),"_",code,sep="")
dbWriteTable(con,result.name,result,row.names=F,overwrite=T)
marks <- unique(sapply(strsplit(colnames(result)[4:length(colnames(result))],"_"),'[',1))
marks <- paste(marks,collapse=";")

if (mode == "Interindividual") {
  sqlquery <- paste("INSERT INTO Interindividual VALUES('"
                    ,result.name,"','",tissue,"','",length(donorlist),"','hg19','Roadmap','",marks,"');",sep="")
  print(tissue,code)
} else if (mode == "Intraindividual") {
  sqlquery <- paste("INSERT INTO Intraindividual VALUES('"
                    ,result.name,"','",donor,"','","Unknown","','",nrow(tissuelist),"','hg19','Roadmap','",marks,"');",sep="")
  print(tissuelist)
}
invisible(dbSendQuery(con,sqlquery))

invisible(dbDisconnect(con))

######################
## FEATURE ANALYSIS ##
######################

# This module of the function attempts to find a correlation between windows with
# epigenetic marks and genomic features such as genes. 

#gff <- read.table('Chr22.gff',header=FALSE)[,c(1,3,4,5)]
#names(gff) <- c("chr","feature","start","end")

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

# https://github.com/rstats-db/RMySQL/issues/140