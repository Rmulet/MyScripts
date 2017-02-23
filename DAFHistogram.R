#!/usr/bin/Rscript

library(ggplot2)

## Total DAF ##

DAF <- t(sapply(strsplit(gen.data$DAF,";"),as.numeric))
# With barplot or similar
freqnames <- c("0.0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1.0")
colnames(DAF) <- freqnames
barplot(colSums(DAF),col="firebrick",border=NA)
barplot(colSums(DAF),col="firebrick",border=NA,space=1,cex.names=0.8)
# With ggplot:
bins <- seq(0.0,1.0,0.1)
d <- data.frame(bins=bins[1:length(bins)-1],DAF=colSums(DAF))
freqnames2 <- c("0.1-0.2","0.3-0.4","0.5-0.6","0.7-0.8","0.9-1.0")
ggplot(data=d,aes(x=bins,y=DAF)) + geom_bar(stat='identity',fill="firebrick") + labs(x="Derived allele frequency",y="Count (in millions of sites)") + scale_x_continuous(breaks=c(0.1,0.3,0.5,0.7,0.9),labels=freqnames2)

# With authentic histogram:

library(HistogramTools) # To create histogram from existing bins
plot(PreBinnedHistogram(bins,colSums(DAF)),col="firebrick",border=NA,ylab="Count (in millions of sites)",xlab="Derived allele frequency",yaxt="n",las=1,main="Unfolded site frequency spectum")
axis(2,at=seq(0,800000,by=200000),label=c(0.0,0.2,0.4,0.6,0.8)) # las = orientation of the labels