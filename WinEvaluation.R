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
  hist(varnum,xlab="Number of variants",breaks=6,main=sprintf("Window %d",wsize),col="steelblue",lty="blank")
  dat.hist <- c(dat.hist,median(varnum))
}
dev.off()
data <- data.frame(Wsize=wsizes,Median=dat.hist)
qplot(y=Median,x=Wsize,data=data)
# plot(dat.hist ~ wsizes,xlab="Window size",ylab="Median number of variants",xaxt="n")
# axis(1,at=wsizes)
qplot(y=Median,x=Wsize,data=data,geom=c("point","smooth"),ylab="Median number of variants",
      main="Number of variants ~ window size") + theme(plot.title = element_text(face="bold",margin = margin(t = 10, b = 10))) +
  theme(axis.title.x = element_text(color="forestgreen",margin=margin(t=10)),axis.title.y = element_text(color="forestgreen",margin=margin(r=10))) + 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) + scale_x_discrete(name ="Window size", limits=wsizes)
