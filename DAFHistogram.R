#!/usr/bin/Rscript

# setwd("D:/MyScripts")
# file <- "DAFTable.tab"

######################
## DATA IMPORTATION ##
######################

# Argument from command line as file
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

# Import DAF file
gen.data <- read.table(file,header=TRUE,stringsAsFactors = FALSE)

# Remove regions with a NA
gen.data <- gen.data[complete.cases(gen.data),]

# Put the DAF in a matrix
DAF <- t(sapply(strsplit(gen.data$DAF,";"),as.numeric))

#####################
## PLOT GENERATION ##
#####################

# Make sure that HistogramTools is installed

if (!"HistogramTools" %in% installed.packages()) {
  cat("HistogramTools was not found on your system and will be installed\n")
  install.packages("HistogramTools",repos="http://cran.r-project.org",INSTALL_opts = c('--no-lock'))
}

library(HistogramTools) # To

# Obtain interval and define bins
interval <- 1/ncol(DAF)
bins <- seq(0.0,1.0,interval)

# Create histogram from existing bins
png('DAFHistogram.png')
plot(PreBinnedHistogram(breaks=bins,counts=colSums(DAF)),col="firebrick",border=NA,ylab="Count (in millions of sites)",xlab="Derived allele frequency",yaxt="n",las=1,main="Unfolded site frequency spectum")
axis(2,at=seq(0,800000,by=200000),label=c(0.0,0.2,0.4,0.6,0.8)) # las = orientation of the labels
dev.off()
  