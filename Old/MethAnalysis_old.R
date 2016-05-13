#!/usr/bin/Rscript

# setwd("~/Documents/3_EpigenomicsData/Roadmap/Methylation")
# file <- "GSM1010978_UCSD.Left_Ventricle.Bisulfite-Seq.STL001.wig"
# file <- "Left_ventricle.BS-Seq.STL001.bw"

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

#######################################
## FUNCTION TO CONVERT WIG TO BIGWIG ##
#######################################  

# WARNING: Requires the use of 'bwtool' and 'WigtoBigWig', downloadable here:
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
# https://github.com/CRG-Barcelona/bwtool/wiki
# For 'WigtoBigWig', chromosome sizes are needed. They can be retrieved with
# the aplication 'fetchChromSizes', available in the first directory.

wigextractor <- function(chr=22,ini=1,step=nrow(win)) { # By default, takes all windows
  library(stringr)
  pattern <- "\\.(\\w+)\\.Bisulfite-Seq.(\\w+)\\."  # Add .wig if used here
  filenames <- list.files(".", pattern="*.wig.gz", full.names=TRUE) # Files in the folder
  samples <- vector()
  for (file in filenames) {
    # EXTRACT INFORMATION FROM THE FILE NAME
    filechr <- paste(substr(file,1,nchar(file)-2),"chr22",sep="")
    out <- paste(substr(file,1,nchar(file)-7),".bw",sep="")
    tissue.id <- str_match(file,pattern)[,2] # Before mark
    donor.id <- str_match(file,pattern)[,3] # After mark
    if (anyNA(c(tissue.id,donor.id))) {
      stop(sprintf("File %s does not follow the required pattern. Please check with '-h' or '--help' \n",file))
    }
    if (mode == "Interindividual" && tissue == tissue.id) {sample <- donor.id}
    else if (mode == "Intraindividual" && donor == donor.id) {sample <- tissue.id}
    else {cat(sprintf("File %s not appropriate for the current analysis \n",file));next}
    samples <- c(samples,sample)
    # IF NECESSARY, CREATE BIGWIG:
    if (file.exists(out) == FALSE) {
      # Select only 'chr22' rows:
      system(sprintf("gunzip -c %s | awk '/chrom=chr%d/{p=1}/chrom=chr[^22]/{p=0}p' > %s",file,chr,filechr))
      # Transform wig to bigWig:
      system(sprintf("wigToBigWig %s hg19.chrom.sizes %s",filechr,out))  
    }
    # EXTRACT INFORMATION IN WINDOWS:
    # Create 'windows.bed' file for 'bwtools':
    write.table(win[ini:(ini+step-1),],file="windows.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    system("sed -i 's/^22/chr22/g' windows.bed") # Only because we don't have "chr"
    # Import as a table to R:
    system(sprintf("bwtool extract bed windows.bed %s %s -decimals=3",out,paste(out,".wn",sep="")))
    tab <- read.table(paste(out,".wn",sep=""),stringsAsFactors = FALSE)
    suppressWarnings(assign(sample,lapply(strsplit(tab[,5],","),as.numeric),envir=.GlobalEnv))
  }  
  return(samples)  
}

############################################
## FUNCTION TO MEASURE WINDOW METHYLATION ##
############################################

methylscore <- function(nwin=nrow(win)) { # By default, it takes all the windows
  # Preallocate vectors containing the results
  meanavg <- numeric(nwin)
  varavg <- numeric(nwin)
  meanvar <- numeric(nwin)
  # Compares the n-line of each window across samples (tissues/donors/species)
  for (window in 1:nwin) {
    mat <- matrix(nrow=0,ncol=200) # Epigenetic diversity matrix
    for (sample in samples) { # We add every sample to the matrix
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
    
######################################
## METHYLATION COMPARATIVE ANALYSIS ##
######################################

if (nrow(win) < 500) {
  samples <- wigextractor()
  methyldata <- methylscore()
} else {
  methyldata <- data.frame(MethLevel=numeric(0),MethLevelVar=numeric(0),MethVar=numeric(0))
  for (i in seq(1,nrow(win),by=100)) {
    samples <- wigextractor(chr=22,ini=i)
    methyldata <- rbind(methyldata,methylscore(nwin=100))
  }
}

if (mode == "Interindividual") {
  code <- abbreviate(str_replace_all(tissue,"_","")) # Summary table
  result.name <- paste(substr(mode,1,7),"_",code,sep="")
  marks <- dbGetQuery(con,sprintf("SELECT EpiMark FROM Interindividual WHERE ExpId='%s'",result.name))
  if (length(marks) != 0) {
  marks <- paste(marks,"Methylation",sep=";")  
  sqlquery <- sprintf("INSERT INTO Interindividual VALUES('%s','%s','Unknown','%s','hg19','Roadmap','%s');",
                      result.name,donor,nrow(tissuelist),marks)
  } else {
    sqlquery <- sprintf("INSERT INTO Interindividual VALUES('%s','%s','Unknown','%s','hg19','Roadmap','Methylation');",
                        result.name,donor,nrow(tissuelist))
  }
} else if (mode == "Intraindividual") {
  code <- donor
  result.name <- paste(substr(mode,1,7),"_",code,sep="")
  marks <- dbGetQuery(con,sprintf("SELECT EpiMark FROM Intraindividual WHERE ExpId='%s'",result.name))
  if (length(marks) != 0) {
    marks <- paste(marks,"Methylation",sep=";")
    sqlquery <- sprintf("UPDATE Intraindividual SET EpiMark='%s' WHERE ExpId='%s'",marks,result.name)
  } else {
    sqlquery <- sprintf("INSERT INTO Intraindividual VALUES('%s','%s','Unknown','%s','hg19','Roadmap','Methylation');",
                        result.name,donor,nrow(tissuelist))
  }
}
invisible(dbSendQuery(con,sqlquery))

if (dbExistsTable(con,result.name) == TRUE) {
  result.tab <- dbReadTable(con,result.name)
  result.tab <- data.frame(result.tab,methyldata)
} else {
  result.tab <- data.frame(windows,methyldata)
}
dbWriteTable(con,result.name,result.tab,overwrite=TRUE,row.names=FALSE)
