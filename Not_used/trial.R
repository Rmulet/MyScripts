filename <- "merge.vcf.gz"; ini <- 31768081-1; end <- 31899603; chrom <- "22"; wsize <- 10000; MK <- TRUE; ref.chr="chr22.fa"
require(PopGenome)
ncolstest <- vector()
for (i in c(100,1000,seq(2000,15000,by=1000))) {
  t1 <- system.time(readVCF(filename,numcols=i,tid=chrom,from=ini,to=end,include.unknown=TRUE,gffpath=sprintf("chr%s.gff",chrom)))[[1]]
  t2 <- system.time(readVCF(filename,numcols=i,tid=chrom,from=ini,to=end,include.unknown=TRUE,gffpath=sprintf("chr%s.gff",chrom)))[[1]]
  t3 <- system.time(readVCF(filename,numcols=i,tid=chrom,from=ini,to=end,include.unknown=TRUE,gffpath=sprintf("chr%s.gff",chrom)))[[1]]
  print(c(i,mean(c(t1,t2,t3))))
  ncolstest <- c(ncolstest,mean(c(t1,t2,t3)))
}