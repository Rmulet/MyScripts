
# setwd("/home/roger/Documents/Genomics/MyScripts")
# # install.packages("PopGenome",repos="http://cran.r-project.org")
# ini=50278879; end=50290610; chrom <- "22"; filename <- "chr22.vcf.gz"

library(PopGenome)

########################
## LOAD GENOME OBJECT ##
########################

GENOME.class <- tryCatch({readVCF(filename,numcols=10000,tid=chrom,from=ini,to=end,include.unknown=TRUE)},error = function(e) {
  message(e); write(sprintf("%s:%d-%d -- %s",chrom,ini,end,e),sprintf("error_chr%s_%s.log",chrom,pop),append=TRUE)
  return(NULL)})

# Verify that the region object is not null (failure to load), contains variants and has been loaded onto R. If it is, then we assume
# that there are no variants, and therefore S,D and all related metrics are 0 (except for alpha)

if (is.null(GENOME.class)||is.logical(GENOME.class)||GENOME.class@n.biallelic.sites==0) { # If readVCF fails, region=FALSE. If no variants, region=NULL
  cat("This region does not contain any variants: DAF set to NA\n\n")
  export <- data.frame(chr=chrom,start=ini,end=end,DAF=NA)
  write.table(export,file=paste(output.name,".tab",sep=""),quote=FALSE,sep="\t",row.names=F,append=TRUE,
              col.names=!file.exists(paste(output.name,".tab",sep="")))  # Column names written if file does not exist
  quit()
}

individuals <- get.individuals(GENOME.class)[[1]]

humans <- individuals[which(individuals != c("Chimp","Chimp.2"))]

n <- length(humans) 

bial <- get.biallelic.matrix(GENOME.class,1) # Biallelic matrix

bialhuman <- bial[1:n,,drop=F] # Remove outgroup (drop = F to keep 1 dimension)

misshuman <- colSums(is.na(bialhuman))>0 # Sites missing in humans (e.g. structural variants)

polym <- apply(bialhuman,2,sum)>0 & !misshuman# Sites polymorphic in humans: non REF (0) alleles

bialhuman <- bialhuman[,polym,drop=F] # Remove sites monomorphic in humans

polyal <- GENOME.class@region.data@polyallelic.sites[[1]] # Positions of all polyalleles

polysites <- length(polyal) # N polyallelic sites in that region

########################
## LOAD GENOME OBJECT ##
########################

m <- GENOME.class@n.sites-sum(misshuman,na.rm=T)-polysites

freqs <- apply(bialhuman,2,table)

k <- sum(freqs[1,]*freqs[2,])

Pi <- function(k,m,n) { 
  
  comb <- choose(n,2) # Binomial coefficient = combination without repetition
  pi <- k/(comb*m)
  
  return(pi)
  
}

Pi(k,m,n)

GENOME.class <- diversity.stats(GENOME.class,pi=TRUE)
GENOME.class@Pi/GENOME.class@n.sites