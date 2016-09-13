##########################
## ALL CODING POSITIONS ##
##########################

CDS <- unique(cbind(START,END,REV))
cds.pos <- apply(CDS,1,function(x){x[1]:x[2]})
all.codons <- lapply(cds.pos,function(x){CHR[x]}) # Select positions specified in codons from FASTA file
synGLOBAL$count <- 1
all.codons <- lapply(all.codons,function(x){
  if (CDS[synGLOBAL$count,3] == 1){
    x <- komplement[x]
  }
  synGLOBAL$count <- synGLOBAL$count+1
  return(x)
})
s <- lapply(all.codons,function(x){
  fold <- vector()
  test <- vector()
  x2 <- x[1:(length(x)-length(x)%%3)] # We remove nt that are not multiple of 3
  indices <- seq(0,length(x2)-1)%%3+1 # +1 because we are using assocmat
  a <- matrix(x2,ncol=3,byrow=TRUE)
  for (i in 1:length(x2)){
    mpos <- ceiling(i/3)
    ncod <- as.numeric(codonise64(a[mpos,,drop=F]))
    test <- c(test,code[ncod])
    fold <- c(fold,assocmat[ncod,indices[i]])
  }
  return(test)
})
