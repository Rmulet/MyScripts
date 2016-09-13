# Loading the R-package parallel
# Lets define the two regions
parallelVCF <- function(name,file,start,end,cores) {
  require(parallel)
  hinge = (end-start)/2
  cregions <- character(2)
  cregions[1] <- paste(start,"-",hinge,sep="")
  cregions[2] <- paste((hinge+1),"-",end,sep="")
  return(cregions)
  
  name.temp <- parallel::mclapply(as.list(cregions),
                             
  function(x){
      From <- as.numeric(strsplit(x,"-")[[1]][1])
      To <- as.numeric(strsplit(x,"-")[[1]][2])
      return(readVCF(
        filename="file",
        numcols=1000, tid="2", frompos=From, topos=To, samplenames=NA,
        gffpath=FALSE, include.unknown=TRUE,
        approx=FALSE, out=x, parallel=FALSE))
  },
  mc.cores = cores, mc.silent = TRUE, mc.preschedule = TRUE)
  
  name <- concatenate.regions(name)
}