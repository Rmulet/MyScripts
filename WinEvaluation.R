#!/usr/bin/Rscript
setwd("~/Documents/2_GenomicsData/TestPopGenome")
filename <- "chr22_filtered.vcf.gz"; ini <- 31768006; end <- 31899628; wsize <- 200
# The centromere of the chr22 is at position 14.7. We avoid it because it has low levels of variation.

suppressMessages(library("PopGenome",quietly=TRUE))

region <- readVCF(filename,numcols=5000,tid="22",from=ini,to=end,include.unknown=TRUE)

individuals <- get.individuals(region)[[1]]
humans <- individuals[which(individuals != c("Chimp","Chimp.2"))]
chimp <- individuals[which(individuals == c("Chimp","Chimp.2"))]
region <- set.populations(region,list(humans,chimp))

wsizes <- seq(from=200,to=2000,by=200)
n <- sqrt(length(wsizes))
par(mfrow=c(n,n+1))
dat.hist <- vector()

for (wsize in wsizes) {
  slide <- sliding.window.transform(region,width=wsize,jump=wsize,type=2)
  varnum <- unlist(lapply(slide@SLIDE.POS,length))
  hist(varnum,xlab="Number of variants",main=sprintf("Window %d",wsize))
  dat.hist <- c(dat.hist,median(varnum))
}
dev.off()
plot(dat.hist ~ wsizes,xlab="Window size",ylab="Median numver of variants",xaxt="n")
axis(1,at=wsizes)
