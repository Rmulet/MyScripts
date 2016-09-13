#!/usr/bin/Rscript

# setwd("~/Documents/3_EpigenomicsData/Roadmap/Methylation")
# file <- "GSM1010978_UCSD.Left_Ventricle.Bisulfite-Seq.STL001.wig"
# file <- "Left_ventricle.BS-Seq.STL001.bw"

################################################
## ACCESS MYSQL DATABASE AND RETRIEVE WINDOWS ##
################################################

headerbed <- c("chr","start","end")

start <- Sys.time()
# WINDOWS IN MYSQL:
suppressMessages(library(RMySQL))
# Retrieve data from MySQL and put in GRanges format
con <- dbConnect(RMySQL::MySQL(),
                 user="root", password="RM333",
                 dbname="PEGH", host="localhost")
res <- dbSendQuery(con, "SELECT Chr,Window FROM Genomics")
windows <- dbFetch(res)
# win <- cbind(windows[,1],trimws(do.call(rbind,strsplit(windows[,2],"-"))))
# Instead of the above, we do this becaused windows are 1-based (PopGenome):
win <- apply(do.call(rbind,strsplit(windows[,2],"-")),2,as.numeric)
win <- cbind(windows[,1],win[,1],win[,2]+1) # Only because PopGenome is 1-based
nwin <- nrow(win)

write.table(win,file="windows.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
system("sed -i 's/^22/chr22/g' windows.bed") # Only because we don't have "chr"

####################################
## MEASURE METHYLATION IN WINDOWS ##
####################################

library(stringr)

pattern <- "\\.(\\w+)\\.Bisulfite-Seq.(\\w+)\\."  # Add .wig if used here
filenames <- list.files(".", pattern="*.wig.gz", full.names=TRUE) # Files in the folder
for (file in filenames) {
  filechr <- paste(substr(file,1,nchar(file)-2),"chr22",sep="")
  out <- paste(substr(file,1,nchar(file)-7),".bw",sep="")
  # Select only 'chr22' rows:
  system(sprintf("gunzip -c %s | awk '/chrom=chr22/{p=1}/chrom=chr[^22]/{p=0}p' > %s",file,filechr))
  # Transform wig to bigWig:
  system(sprintf("wigToBigWig %s hg19.chrom.sizes %s",filechr,out))
  # Extract windows from the bigWig file:
  system(sprintf("bwtool extract bed windows.bed %s %s -decimals=3",out,paste(out,".wn",sep="")))
}

filenames2 <- list.files(".", pattern="*.wn", full.names=TRUE) # Files in the folder
samples <- vector()
for (file in filenames2) {
  tissue.id <- str_match(file,pattern)[,2]
  donor.id <- str_match(file,pattern)[,3]
  if (anyNA(c(tissue.id,donor.id))) {
    stop(sprintf("File %s does not follow the required pattern. Please check with '-h' or '--help'",file))
  }
  tab <- read.table(file,stringsAsFactors = FALSE)
  suppressWarnings(assign(tissue.id,lapply(strsplit(tab[,5],","),as.numeric)))
  samples <- c(samples,tissue.id)
}
start <- Sys.time()
# Preallocate vectors containing the results
meanavg <- numeric(nwin)
varavg <- numeric(nwin)
meanvar <- numeric(nwin)
for (window in 1:nwin) {
  mat <- matrix(nrow=0,ncol=200) # Epigenetic diversity matrix
  for (sample in samples) { # We add every sample to the matrix
    mat <- rbind(mat,eval(as.symbol(sample))[[window]])
  }
  matmet <- mat[,complete.cases(t(mat))] # We remove NA (non-meth/missing)
  avg <- apply(matmet,1,mean)
  meanavg[window] <- mean(avg)
  varavg[window] <- var(avg)
  posvar <- apply(matmet,2,var)
  meanvar[window] <- mean(posvar)
}
methyldata <- data.frame(MethLevel=meanavg,MethLevelVar=varavg,MethVar=meanvar)

Sys.time() - start

