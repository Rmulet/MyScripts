#!/usr/bin/Rscript

# setwd("~/Documents/2_GenomicsData/Final/GeneByGene")
# gpfile <- "chr21_gp.vcf.gz" ; alnfile <- "chr21_aln.vcf.gz"; chrom <- "21" ; maskfile <- "20140520.chr21.pilot_mask.fasta.gz"

# UPDATE: Instead of segmenting the PopGenome object, we will remove those positions that
# are missing via vector.
# UPDATE 2: PopGenome steps have been separated as a distinct function.
# UPDATE 3: It is no longer necessary to import FASTA or GFF, as we are using ourw own vector.
# UPDATE 4: Memory limitations due to excessive size of bialhuman and bial.
# UPDATE 5: Correction in the way m/m0 are calculated. 
# UPDATE 6: File naming updated to allow parallelization.

suppressMessages(library(PopGenome))
suppressMessages(library(GenomicRanges))
suppressMessages(library(stringr))
suppressMessages(library(DBI))
suppressMessages(library(RMySQL))

# zz <- file("Warnings.dat","w")
# sink("Warnings.dat",append=TRUE)

#############################
## IMPORT AND PREPARE DATA ##
#############################

args <- commandArgs(trailingOnly = TRUE) # Import arguments from command line
gpfile <- args[1] # 1000GP data
alnfile <- args[2] # Chimp-Human alignment
maskfile <- args[3] 
chrom <- args[4]

if (args[1] == "-h" | args[1] == "--help") {
  cat("\nGFFtoFASTA - A script that creates a FASTA-like file with numeric codes identifying each feature.\n")
  cat("\nUsage: GFFtoFASTA.R [GP FILE] [ALN FILE] [MASKFILE]\n\n")
  quit()
}

nfields <- 87 # Number of columns: 7 for polymorphism + 5*16 for selection

##########################
## FUNCTION DECLARATION ##
##########################

## POLYMORPHISM FUNCTIONS ##

Pi <- function(k,m,n) { # Window number
  comb <- choose(n,2) # Binomial coefficient = combination without repetition
  pi <- k/(comb*m)
  return(round(pi,7))
}
Theta <- function(S,m,n) { # S=Segregating sites; m=total sites; n=number of samples
  summ <- 0
  for (i in 1:(n*2-1)) {
    summ <- summ + 1/i
  }
  theta <- (S/m)/summ
  return(round(theta,7))
}

## MKT CALCULATION ##

mkt.extended <- function (sel=0,neu=4,gffseq) {
  
  # Divtotal, poly.sites and bial.class are in Global to pass them
  
  # STANDARD: Same result as MKT function if done with all syn.
  Pneu <- sum(poly.sites[bial.class==neu],na.rm=TRUE) # Psel + Dsel < total sel because of sites that are polym and divtotal
  Psel <- sum(poly.sites[bial.class==sel],na.rm=TRUE)
  Dneu <- sum(divtotal[bial.class==neu],na.rm=TRUE)
  Dsel <- sum(divtotal[bial.class==sel],na.rm=TRUE)
  NI <- (Psel/Pneu)/(Dsel/Dneu)
  alpha <- 1-NI
  DoS <- Dsel/(Dsel+Dneu)-Psel/(Psel+Pneu) # Direction of selection
  contingency.std <- matrix(c(Psel,Pneu,Dsel,Dneu),c(2,2))
  test <- if(!is.na(sum(contingency.std))){
    fisher.test(contingency.std)$p.value
  } else {NA}
  
  # ADAPTED:
  
  # 1) Use four-fold degenerated sites 2) Remove sites according to DAF
  
  neuMAF <- as.vector(na.omit(MAF[bial.class==neu & poly.sites==TRUE]))
  selMAF <- as.vector(na.omit(MAF[bial.class==sel & poly.sites==TRUE]))
  
  Pneu.less5 <- sum(neuMAF<0.05) # P0 MAF less than 5%
  Pneu.more5 <- sum(neuMAF>0.05) # P0 MAF more than 5% 
  Psel.less5 <- sum(selMAF<0.05) # Pi MAF less than 5% 
  Psel.more5 <- sum(selMAF>0.05) # Pi MAF more than 5%
  
  Psel.neutral.less5 <- Psel*(Pneu.less5/Pneu) # Proportion of neutral within the MAF < 5% class
  Psel.neutral <- Psel.neutral.less5 + Psel.more5 # For alpha
  Psel.weak <- Psel.less5 - Psel.neutral.less5
  
  alpha.cor <- 1-(Psel.neutral/Pneu)*(Dneu/Dsel)
  contingency <- matrix(c(Psel.neutral,Pneu,Dsel,Dneu),c(2,2))
  test.cor <- if(!is.na(sum(contingency))){
    if(sum(contingency)>0){fisher.test(contingency)$p.value}
  } else {NA}
  
  # OTHER ESTIMATORS:
  
  # To fully implement extended MKT, we need ms/mns, that is, the number of sites
  # of each class. Doing that would require modifying the 'set.synnonsyn2' function
  # to obtain all codons and check fold in every position (0,1,2)
  
  m.neu <- sum(gffseq[ac.pos] == neu,na.rm=T)
  m.sel <- sum(gffseq[ac.pos] == sel,na.rm=T)
  
  f <- (m.neu*Psel.neutral)/(m.sel*Pneu) # Neutral sites
  b <- (Psel.weak/Pneu)*(m.neu/m.sel)
  y <- (Psel/Pneu-Dsel/Dneu)*(m.neu/m.sel)
  d <- 1 - (f+b)
  
  return(c(Psel,Pneu,Dsel,Dneu,m.neu,m.sel,alpha,test,Psel.neutral,alpha.cor,test.cor,DoS,f,b,y,d))
}

#################################
## POPGENOME ANALYSIS FUNCTION ##  
#################################

popanalysis <- function(filename,ini,end,chrom,ac.pos,gffseq) {
  vcount <- as.numeric(system(sprintf('zcat %s | grep -v "#" | wc -l',filename),intern=TRUE))
  if (vcount == 0) { # Check if there are variants to prevent PopGenome error
    cat ("There are no variants in this VCF file")
    newrow <- c(rep(0,6),rep(NA,81)) # Empty rows
    return(newrow)
  }
  region <- readVCF(filename,numcols=9000,tid=chrom,from=ini,to=end,include.unknown=TRUE)
  # Syn-nonsyn is not needed if we only use 0- and 4-fold. Therefore, GFF and FASTA can be skipped.

  ## DEFINE POPULATION/OUTGROUP ##
  
  # The outgroup can be defined as a different population so that calculations of diversity and neutrality 
  # are performed separately. 
  
  individuals <- get.individuals(region)[[1]]
  humans <- individuals[which(individuals != c("Chimp","Chimp.2"))]
  chimp <- individuals[which(individuals == c("Chimp","Chimp.2"))]
  region <- set.populations(region,list(humans,chimp))
  region <- set.outgroup(region,chimp)
  
  ## FILTER THE BIALLELIC MATRIX ##
  
  n <- length(humans) # N = samples. Diploid, not divided by 2
  total <- length(individuals) # Total number of samples (including outgroup)
  bial <- get.biallelic.matrix(region,1) # Biallelic matrix of the window
  ac.bial <- as.numeric(colnames(bial)) %in% ac.pos # Accessible positions in biallelic matrix
  bial <- bial[,ac.bial,drop=F] # Biallelic matrix of accessible variants
  if (is.null(bial)||dim(bial)[2]==0) { # When no variants are detected
    newrow <- c(rep(0,6),rep(NA,81)) # Empty rows
    return(newrow)
  }
  wsize <- end-ini+1 # GFF coordinates, 1-based
  
  # Bial contains only positions that are found in the accessibility mask
  
  misshuman <- colSums(is.na(bial[1:n,,drop=F]))>0 # Sites missing in humans (e.g. structural variants)
  poly.sites <<- apply(bial[1:n,,drop=F],2,sum)>0 & !misshuman # Sites polymorphic in humans w/o missing
  #bialhuman <- bialhuman[,poly.sites,drop=FALSE] # Keep only polymorphic sites
  #bial[1:n,poly.sites,drop=F] # Remove outgroup (drop = F to keep 1 dimension)
  
  cat("\nBiallelic matrix filtered\n")
  
  ## SITE FREQUENCY SPECTRUM (FOLDED) ##
  
  MAF <<- colSums(bial[1:n,,drop=F] == 1)/n # We DO NOT remove monomorphic alleles
  MAF.df <- data.frame(POS=names(MAF),MAF=unname(MAF)) # Position and DAF of biallelic variants
  
  # IDENTIFY ANCESTRAL ALLELES (AA) FOR SNPS IN THE REGION
  temp <- system(sprintf('gunzip -c %s | cut -f2,4,8',filename),intern=TRUE) # Extract all variants from the merge.vcf.gz file  
  gpimport <- read.table(textConnection(temp),sep="\t",stringsAsFactors = FALSE,header=TRUE)
  gpances <- toupper(sapply(gpimport[,3],function(x){str_match(x,"AA=([:alpha:])\\|")[2]},USE.NAMES=FALSE)) # Capture the ancestral allele
  aarefs <- na.omit(cbind(gpimport[,1:2],ANC=gpances)) # Na.omit removes positions without known AA and divergent sites
  
  # CALCULATE DERIVED ALLELE FREQUENCY (DAF)
  mafan <- merge(MAF.df,aarefs) # Intersection between SNPs in the biallelic matrix and the AA in the VCF file
  if (dim(mafan)[1] > 0) { # Avoid windows where variants have no known AA
    mafan$DAF <- NA # DAF is equal to MAF if derived = minor, but 1 - MAF if otherwise
    mafan[mafan$REF != mafan$ANC,]$DAF <- 1- mafan[mafan$REF != mafan$ANC,]$MAF 
    mafan[mafan$REF == mafan$ANC,]$DAF <- mafan[mafan$REF == mafan$ANC,]$MAF 
  }
  mafan <- merge(MAF.df,mafan,all.x=TRUE) # Retain originals
  
  if (is.numeric (mafan$DAF)) {
    DAF <- hist(mafan$DAF,seq(0.0,1,0.05),plot=F)$counts
    DAF <- paste(DAF,collapse = ";")
  } else {DAF <- NA}
  Sys.sleep(3)
  
  print(sort(sapply(ls(),function(x){object.size(get(x))}))) # Memory of each object
  
  ## POLYMORPHISM ## 
  
  # DETERMINE M (EXCLUDING MISSING AND POLYALLELIC)
  # Total number of sites: 1) Remove NA in humans 2) Remove polyallelic sites
  if (!length(region@region.data@polyallelic.sites) == 0) { # Make sure list exists
    multisites <- length(region@region.data@polyallelic.sites[[1]]) # Positions of all polyalleles
     } else {multisites <- 0} # If not available, we assume 0
  m <- wsize-sum(misshuman,na.rm=T)-multisites
  # DETERMINE S AND K (WHEN VARIANTS ARE AVAILABLE)
  if (is.null(bial[1:n,poly.sites,drop=F])||dim(bial[1:n,poly.sites,drop=F])[2]==0) {
    S <- 0
    k <- 0
  } else {
    S <- ncol(bial[1:n,poly.sites,drop=F]) # Number of variants
    freqs <- apply(bial[1:n,poly.sites,drop=F],2,table)
    k <- sum(freqs[1,]*freqs[2,]) # Note that k and K are different!
  }
  
  # Functions are declared above!
  
  ## DIVERGENCE ##
  
  # COUNT THE NUMBER OF UNKNOWNS IN THE OUTGROUP:
  unknowns <- sum(is.nan(bial[total,])) # Displayed as NaN in the biallelic matrix
  
  # COUNT THE NUMBER OF DIVERGENCE SITES:
  divtotal <<- bial[total,] == 1 & !poly.sites # 1/1 in the outgroup and excluding polymorphic (cannot be counted for divergence)
  divsites <- sum (divtotal,na.rm=T) # Number of divergent sites
  # CALCULATE D AND K:
  mout <- m - unknowns
  D <- round(divsites/mout,7) # Observed divergence: Proportion of sites with divergent nucleotides
  K <- round(-3/4*log(1-4/3*D),7) # Real divergence: Jukes and Cantor model to account for multiple hits
  
  cat("Divergence and polymorphism calculated\n")
  
  ## NATURAL SELECTION REGIMES ##
  
  bial.sites <- as.numeric(colnames(bial)) # Biallelic sites
  bial.class <<- gffseq[bial.sites] # GFF feature in each site. Equivalent to the former 'syn'
  # syn <- region@region.data@synonymous[[1]][ac.bial] # Contains syn and non-syn positions [accessible]
  
  # STORING MKT RESULTS FOR DIFFERENT FUNCTIONAL CLASSES:
  
  mkt.0fold.4fold <- mkt.extended(sel=0,neu=4,gffseq) # sel = 0-fold; neu= 4-fold
  mkt.intron.4fold <-  mkt.extended(sel=9,neu=4,gffseq) # sel = exon; neu= 4-fold
  mkt.5UTR.4fold <-  mkt.extended(sel=5,neu=4,gffseq) # sel = 5-UTR; neu= 4-fold
  mkt.3UTR.4fold <-  mkt.extended(sel=3,neu=4,gffseq) # sel = 3-UTR; neu= 4-fold
  bial.class[is.na(bial.class)] <- 1 # New code for all intergenic (instead of NA)
  mkt.inter.4fold <-  mkt.extended(sel=1,neu=4,gffseq) # sel = exon; neu= 4-fold
  
  cat("MKT performed. Data will be stored in a new row\n")
  
  ## ADD NEW ROW ##
  newrow <- c(S,Pi(k,m,n),DAF,divsites,D,K,unknowns,mkt.0fold.4fold,mkt.intron.4fold,mkt.5UTR.4fold,mkt.3UTR.4fold,mkt.inter.4fold)
  memory.profile()
  return(newrow)  
}

###################
## GENE ANALYSIS ##
###################

merge.vcf <- function(ini,end,filename) { # Indicate the directory if experimental bcftools is local (Andromeda)
  t <- try(system(sprintf("bcftools merge -Oz --missing-to-ref -o %s -r %s:%d-%d %s %s",
                          filename,chrom,ini,end,gpfile,alnfile)))
  if ("try-error" %in% class(t)) {
    gc(reset=T)
    system(sprintf("bcftools/bcftools merge -Oz --missing-to-ref -o %s -r %s:%d-%d %s %s",
                   filename,chrom,ini,end,gpfile,alnfile))
  }
  system(sprintf("tabix -p vcf %s",filename))
}

## RETRIEVE ENTREZ GENES AND PREPARE DATA TABLE:

load("GenesTable.RData")
gendata <- genestable[genestable$chr == sprintf("chr%s",chrom),] # Select genes in this chromosome
gendata <- gendata[order(gendata$start),]
gendata$start <- gendata$start-500 # We chose to sequence 500 bp upstream and downstream [PMID: 21059791
gendata$end <- gendata$end+500
ngenes <- nrow(gendata)

# TABLE WITH DATA:
tabsum <- as.data.frame(matrix(numeric(ngenes*nfields),ncol=nfields,nrow=ngenes))
colnames <- paste(c("S","Pi","DAF","Divsites","D","K","Unknown"),rep(c("Psel","Pneu","Dsel","Dneu","alpha","test","Psel.neutral","alpha.cor","test.cor","DoS","f","b","y","d"),6))

## PREANALYSIS: GFF AND MASK:

load(sprintf("gffseq_chr%s.RData",chrom))

library("Biostrings")
maskfasta <- readBStringSet(maskfile) # Reading files in gz format IS supported
cat("Maskfile loaded\n")

# Note that the coordinates of the MASK are in BED format and therefore 0-based, whereas the GFF with
# the genes and the GRanges objects are 1-based. To convert from 0 to 1-based, START+1:END.

## GENE BY GENE ANALYSIS:

init <- Sys.time()
for (i in 1:20) {
  ini <- gendata$start[i]; end <- gendata$end[i]
  print(sprintf("Gene number %d: %d - %d",i,ini,end))
  mask.local <- strsplit(as.character(subseq(maskfasta,start=ini,end=end)),"")[[1]]
  pass <- mask.local == "P"
  print(sum(pass))
  if (sum(pass) == 0) {
    gendata$missing[i] <- 100
    tabsum[i,] <- rep(NA,nfields) 
    next}
  ac.pos <- (ini:end)[pass] # Vector with gene positions that are accessible
  gendata$missing[i] <- (1-sum(pass)/length(pass))*100 # Proportion of positions that do not                                                                                                                                                                                                                   pass the filter
  filename <- sprintf("merge_gene%s",gendata$name[i])
  merge.vcf(ini,end,filename)
  tabsum[i,] <- popanalysis(filename,ini,end,chrom,ac.pos,gffseq)
  # print(sort(sapply(ls(),function(x){object.size(get(x))})))
  if(end-ini > 500000) {gc(reset=T)}
  file.remove(filename)
}
Sys.time()-init

###################
## FINAL TOUCHES ##
###################

# tabsum[,-3] <- apply(tabsum[,-3],2,as.numeric) # Convert all but DAF to numeric
db.names <- c("S","Pi","DAF","Divsites","D","K","Unknown")
mkt.names <- c("Psel","Pneu","Dsel","Dneu","alpha","test","Psel_neutral","alpha_cor","test_cor","DoS","f","b","y","d")
site.class <- c("fold0","intron","UTR5","UTR3","UTR","inter")

for (class in site.class) {
  db.names <- c(db.names,unname(sapply(mkt.names,function(x){paste(c(x,"_",class),collapse="")})))
}
colnames(tabsum) <- db.names

#library(biomaRt)
#mart <- useMart(host="www.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL")
#atts <- c("entrezgene","hgnc_symbol") # Retrieve the HUGO-approved gene symbol
#gene.names <- genes$gene_id
#annotation.h <- getBM(attributes=atts,uniqueRows=FALSE,filters="entrezgene",values=gene.names,mart=mart)

######################
## DATA EXPORTATION ##
######################

# We add the chromosome info and the coordinates. For convenience in the next steps:
# a) Coordinates are converted to BED format: 0-based, subtract 1 from start column.
# b) Chromosome is expressed in "chrNN" format.

export <- cbind(gendata,tabsum)
export.name <- "GenesDB"

# EXPORT DATA TABLE
write.table(export,file=sprintf("GeneData_chr%s.tab",chrom),quote=FALSE,sep="\t",row.names=F)

# EXPORT TO MYSQL DATABASE
con <- dbConnect(RMySQL::MySQL(),
                 user="root", password="RM333",
                 dbname="PEGH", host="localhost")
first <- !dbExistsTable(con,export.name)

t <- try(dbWriteTable(con,value=export,name=export.name,row.names=F,append=T))

if ("try-error" %in% class(t)) {
  save(export,file=sprintf("failedexport_chr%s.RData",chrom))
}

if (first == TRUE) {# Remove if we want to concatenate various chromosomes
  dbSendQuery(con,sprintf("ALTER TABLE %s CHANGE COLUMN name name VARCHAR(30);",export.name))
  dbSendQuery(con,sprintf("ALTER TABLE %s ADD PRIMARY KEY (name);",export.name))
}

on.exit(dbDisconnect(con))
closeAllConnections()
