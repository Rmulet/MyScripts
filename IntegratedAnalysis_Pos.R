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
res <- dbSendQuery(con, "SELECT chr,start,end FROM Genomics")
windows <- dbFetch(res)
invisible(dbClearResult(res)) # Frees resources associated with the query
windows[,2:3] <- apply(windows[,2:3],2,as.numeric) # Make the variable numeric
wsize <- windows[1,3]-windows[1,2] # Windows size
chr <- substr(windows[1,1],4,nchar(windows[1,1])) # Chromosome number
ntotal <- nrow(windows) # Number of windows

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
file.remove(filechr)
return(out)
}

######################################
## FUNCTION TO EXTRACT DATA FROM BW ##
######################################  

# EXTRACT INFORMATION IN WINDOWS:
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
    # Import as a table to R:
    system(sprintf("bwtool extract bed windows.bed %s %s -decimals=3",file,paste(file,".wn",sep="")))
    tab <- read.table(paste(file,".wn",sep=""),stringsAsFactors = FALSE)
    # WARNING: The table is imported to a variable called SAMPLE.
    suppressWarnings(assign(sample,lapply(strsplit(tab[,5],","),as.numeric),envir=.GlobalEnv))
    samples <- c(samples,sample)
  }
  return(samples)
}
# WE CAN INTEGRATE THIS IF THERE IS ONLY ONE SEGMENT OF WINDOWS   

################################################
## FUNCTION TO MEASURE EPIGENETIC VARIABILITY ##
################################################

episcore <- function(nwin=ntotal) { # By default, it takes all the windows
  # Preallocate vectors containing the results
  meanavg <- numeric(nwin)
  varavg <- numeric(nwin)
  meanvar <- numeric(nwin)
  # Compares the n-line of each window across samples (tissues/donors/species)
  for (window in 1:nwin) {
    mat <- matrix(nrow=0,ncol=wsize) # Epigenetic diversity matrix
    for (sample in samples) { # We add every sample to the matrix
      mat <- rbind(mat,eval(as.symbol(sample))[[window]])
    }
    matclean <- mat[,complete.cases(t(mat))] # We remove positions with NA (non-meth/missing)
    if (dim(matclean)[2] == 0) { # No methylated positions in this region
      meanavg[window] <- 0
      varavg[window] <- 0
      meanvar[window] <- 0
      next
    }
    avg <- apply(matclean,1,mean) # Average of methylation in METHYLATION LOCI
    meanavg[window] <- mean(avg)
    varavg[window] <- var(avg)
    posvar <- apply(matclean,2,var)
    meanvar[window] <- mean(posvar)
  }
  epidata <- data.frame(Level=meanavg,LevelVar=varavg,Var=meanvar)
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
    samples <- bwfraction()
    markdata <- episcore()
  } else { # If the region contains > 500 kb, the analysis is fractioned
    markdata <- data.frame(Level=numeric(0),LevelVar=numeric(0),Var=numeric(0))
    chunk <- 500000/wsize
    for (i in seq(1,ntotal,by=chunk)) { # REVIEW INTERVAL!!!!
      if ((n-i+1) < chunk) {
        # Create 'windows.bed' file for 'bwtools':
        chunk <- ntotal%%chunk
      }
      bwfraction(ini=i,step=chunk)
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
 