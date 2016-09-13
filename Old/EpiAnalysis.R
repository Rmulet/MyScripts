setwd("~/Documents/3_EpigenomicsData/Roadmap")
file <- "a.txt"
# grep "chr22" BI.Adipose_Nuclei.H3K4me1.7.narrowPeak | cut -f1,2,3 |  sort > a.txt

source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
library("GenomicRanges")

headerbed <- c("chr","start","end")
headerpeaks <- c("chr","start","end","name","score","strand","signalValue","pValue","qValue","peak")

peaks <- read.table(file,header=F,stringsAsFactors=F)
names(peaks) <- headerbed
grpeaks <- with(peaks,GRanges(chr,IRanges(start,end)))
require(RMySQL)
# Retrieve data from MySQL and put in GRanges format
con <- dbConnect(RMySQL::MySQL(),
                 user="root", password="RM333",
                 dbname="PEGH", host="localhost")
res <- dbSendQuery(con, "SELECT Chr,Window FROM Genomics")
windows <- dbFetch(res)
dbClearResult(res)
coords <- trimws(do.call(rbind,strsplit(windows[,2],"-"))) # Separate the coordinates bound by "-"
coords <- apply(coords,2,as.numeric) # Make the variable numeric
windows <- data.frame(paste("chr",windows[,1],sep=""),coords)
names(windows) <- headerbed
grwindows <- with(windows,GRanges(chr,IRanges(start,end)))

# FIND OVERLAPS BETWEEN WINDOWS AND CHIP-SEQ DATA:
ov.index <- findOverlaps(query=grwindows,subject=grpeaks)
span <- width(ranges(ov.index,ranges(grwindows),ranges(grpeaks))) # Measures span of the overlap
values <- span/200*100 

query <- queryHits(ov.index)
dup <- query[duplicated(query) | duplicated(query,fromLast=TRUE)]
duppos <- which(duplicated(query) | duplicated(query,fromLast=TRUE))
duppeaks <- subjectHits(ov.index)[duppos]

a <- c(400,402,402,402,403,404,405,405,406) # Equivalent to query
which <- which(duplicated(a) | duplicated(a,fromLast=TRUE))
# dpeakies <- c(1847,1848,1849,1890,1891)
# adup <- a[which]
adup <- unique(a[duplicated(a)])
val <- c(100,20,30,20,67,23,12,31,69)
totalscores <- vector()
for (i in adup) {
  score <- 0
  rep <- which(a == i)
  valrep <- val[rep]
  peakspos <- subjectHits(ov.index)[rep]
  tab <- peaks[peakspos,]
  tabwin <- windows[i,]
  for (n in 2:length(rep)) {
      if (tab[n,2] > max(tab[1:n-1,3])) {
        score <- valrep[n-1]+valrep[n]
      }
      if (tab[n,2] < max(tab[1:n-1,3])) {  
        if (tab[n,3] < max(tab[1:n-1,3])) {
          next
        }
        else if (tab[n,3] > max(tab[1:n-1,3])) {
          score <- valrep[n-1]+(valrep[n]-(tab[n,3]-max(tab[1:n-1,3])))
        }
      }
  }
  totalscores <- c(totalscores,score)    
}  
dbDisconnect(con)