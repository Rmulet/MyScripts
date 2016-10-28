#!/bin/bash

# GENOME ANALYZER 0.2 - Calls the VCFmerger.sh script for the entire genome

# UPDATE 20161025: Feature for specific chromosome analysis added

# Add DOWNLOAD and MASK option

display_usage() { 
	echo -e "\nThis script analyses patterns of variation along the entire genome by repeatedly calling VCFmerger.sh for each chromosome" 
	echo -e "\nUsage:\n $(basename "$0") [-w --window] [-msk --mask] [-dl --download]"
	echo -e "\nOptions:\n -w, --window \t Specifies the size of the window to be analyzed [10000]" 
	echo -e " -msk, --mask \t Indicates what 1000 Genomes Project mask should be used (Pilot/Strict) [Pilot]"
	echo -e " -dl, --download \t Determines whether the necessary files should be downloaded from predefined URLs [FALSE]"
	echo -e " -pop, --population \t Conducts the analysis by population (26/5/name/FALSE) [FALSE]"
	echo -e " -chr, --chromosome \t Analyses a single chromosome (in format NN, e.g. 22) [ALL]"
	} 

#####################################
## ARGUMENT EVALUATION AND PARSING ##
#####################################

## ARGUMENT EVALUATION ##

# Check whether user had supplied -h or --help . If yes display usage 
if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
then 
	display_usage
	exit 0
fi 	

## ARGUMENT PARSING ##

WINDOW=10000
MASK="pilot"
POP="FALSE"
DL="FALSE"
CHR=`seq 1 22` # ALL chromosomes

while [[ $# -gt 1 ]]
do
	case "$1" in
		-w|--window)
		WINDOW="$2" # $1 has the name, $2 the value
		echo -e "Window size set to $WINDOW"
		shift # next two arguments (window + size)
		;;
		-msk|--mask)
		MASK=$(echo "$2" | tr '[:upper:]' '[:lower:]') # $1 has the name, $2 the value
		echo -e "1000 GP mask set to $MASK"
		shift
		;;
		-dl|--download)
		DL="$2" # $1 has the name, $2 the value
		echo -e "ALL necessary files will be downloaded: $DL"
		shift
		;;
		-pop|--population)
		POP="$2" # $1 has the name, $2 the value
		if [ "$POP" == "26" ]; then echo -e "All populations will be analysed\n"
		elif [ "$POP" == "5" ]; then echo -e "The 5 super-populations will be analysed\n"
		else echo -e "Population $POP will be analysed\n"; fi
		shift
		;;		
		-chr|--chromosome)
		CHR="$2" # $1 has the name, $2 the value
		echo -e "Only chromosome $CHR will be analysed\n"
		shift
		;;
		*) # No more options
		;;
	esac
shift
done

##############################
## VARIABLES AND DATA PATHS ##
##############################

WORKING="~/Genomics"

gpraw="$WORKING/1000GP/Chromosomes" # VCF files from 1000 GP divided by chromosomes
gpdat="$WORKING/1000GP" # No files required 
alnraw="$WORKING/Alns/Chromosomes" # Human-chimp alignment (MFA.GZ) divided by chromosomes
alndat="$WORKING/Alns" # Contains FASTA files (FA.GZ/FA)
finaldir="$WORKING/Final" # Contains GFF files

mskfile=$gpdat/Masks/$(cd $gpdat/Masks/ && ls -d *$MASK\_mask.whole_genome.bed) # Depends on the chosen criteria. Underscore must be escaped.

if [ "$DL" == "TRUE" ]; then # DOWNLOAD?
	echo "The files required for the analysis will be downloaded"	
	cd $gpraw
	wget -nc -nd -r -l 1 -A "ALL.chr*" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
	cd $gpdat/Others
	wget -nc ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel -O "PopulationIndividualsList.panel"
	cd $alnraw
	wget -e robots=off -nc -nd -r -l1 -np -A chr$i.mfa.gz,chr$i.mfa.gz http://pipeline.lbl.gov/data/hg19_panTro4/ # Alignments - VISTA Browser 
	wget -e robots=off -nc -nd -r -l1 -np -A chr$i.fa.gz,chr$i.fa.gz ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/
	cd $gpdata/Masks
	wget -nd -r -l0 -np -A "*mask.whole*" -R "index.html*","*combined*" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/
fi	

genome_analysis() {
	for i in $CHR; do
		# POLYMORPHISM #
		echo -e "Processing polymorphism data: chr$i" 
		gpfile=chr$i\_gp.vcf.gz # Name of the filtered file
		cd $gpdat
		if [ ! -e chr$i\_gp.vcf.gz ]; then # Only if the files have not been previously filtered. Not quotes because it's not a variable.
			echo -e "Filtering out inbred individuals in chr$i"
			cd $gpraw
			gpstart=$(ls -d *chr$i.*vcf.gz) # Starting file
			bcftools view -Oz --force-samples -S ^InbredIndividuals.txt $gpstart > $gpfile # Filter out inbred individuals
			mv $gpfile $gpdat; cd $gpdat
			tabix -p vcf $gpfile
		fi

		if [ "$i" == "X" ]; then # Remove MALES from the X chromosome file
			echo -e "Excluding males from chromosome X"
			fem=$(cd $gpdat/Others && grep "female" PopulationIndividualsList.panel | cut -f1 | tr '\n' ',')
			bcftools view -Oz --force-samples -s $fem $gpfile > chr$i.temp.vcf.gz # Remove female individuals
			mv chr$i.temp.vcf.gz $gpfile
			tabix -p vcf chr$i.temp.vcf f $gpfile
		elif [ "$i" == "Y" ]; then
			echo -e "Converting chromosome Y to pseudo-diploid"
			zcat $gpfile | perl -p -e 's/\t([01.])(?=[\n|\t])/\t\1\|\1/g' > chr$i.temp.vcf # Duplicate haploid individuals
			#echo '##fileformat=VCFv4.1\n##FILTER=<ID=PASS,Description="All filters passed">"' > chr$i.temp.vcf
			#bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\tAA=%INFO/AA;AC=%INFO/AC;AN=%INFO/AN\tGT\t[\t%GT|%GT]\n' >> chr$i.temp.vcf
			bgzip chr$i.temp.vcf; mv chr$i.temp.vcf.gz $gpfile
			tabix -p vcf $gpfile
		fi	

		if [ "$POP" != "FALSE" ]; then # Note that "PopulationIndividualsList.panel" is assumed to integrated_call_samples_v3ontain 
			popname=$1
			echo "1"
			popins=$(cd $gpdat/Others && grep $popname PopulationIndividualsList.panel | cut -f1 | tr '\n' ',' | sed 's/,$//' )  # List of individuals in that population
			echo "Extracting the individuals of the selected population: $popname"
			bcftools view -Oz --force-samples -s $popins chr$i\_gp.vcf.gz > chr$i$1.vcf.gz # Some have been removed because they are inbred (force-samples to skip)
			gpfile=chr$i$1.vcf.gz # In population mode, then gpfile is the 
			echo $gpfile	
			tabix -p vcf $gpfile
		fi

		 # DIVERGENCE #
		echo -e "Processing divergence data: chr$i" 
		alnfile=chr$i\_aln.vcf.gz # Alignment file in VCF
		cd $alndat
		if [ ! -e "$alnfile" ]; then
			cd $alnraw
			alnstart=chr$i\_aln.mfa # Alignment file in MFA	
			gunzip -c chr$i.mfa.gz | grep -v "score" > $alnstart # Remove 'score' lines from the MFA file
			mv $alnstart $alndat; cd $alndat
			gunzip chr$i.fa.gz # Uncompress the FASTA sequence
			echo -e "Converting .MFA alignment to VCF for chr$i"
			MFAtoVCF.py -s chr$i.fa -q Chimp -c $i $alnstart # Convert MFA to VCF. VCF.GZ file is automatically tabixed.
		fi

		# PREANALYSIS (GFF to FASTA) #	
		cd $finaldir
		echo -e "Extracting the annotation file to a sequence format"
		cp $alndat/chr$i.fa $finaldir # Copy the FASTA sequence of the chromosome
		if [ ! -e "gffseq_chr$i.RData" ]; then
		GFFtoFASTA8.R chr$i.fa $i # Puts the GFF annotation in a sequence
		fi			

		# MERGE AND ANALYSIS #
		echo -e "Preparing the accessibility mask for chr$i"
		grep "chr$i" $mskfile > $MASK\_mask.chr$i.bed # Extract the chromosome of interest from the mask.		
		ln -s $gpdat/$gpfile $gpfile;  ln -s $alndat/$alnfile $alnfile # Create symbolic links for the 1000GP and Human-Chimp data
		ln -s $gpdat/$gpfile.tbi $gpfile.tbi;  ln -s $alndat/$alnfile.tbi $alnfile.tbi # Create symbolic links for the index files
		echo -e "Analysing polymorphism and divergence in chr$i"
		VCFmerger.sh $MASK\_mask.chr$i.bed $gpfile $alnfile -w $WINDOW -db Genomics$i$popname # Merges and analyzes variation data

		if [[ $? -ne 0 ]]; then # Stop the execution of the script if VCFmerger fails
			echo 'Error: VCFmerger.sh failed!'
			exit -1
		fi

		# The $popname variable will be ommitted if it is not declared (i.e. analysis conducted for 2,504 individuals)
		# By default, the results will be stored in separate tables. To merge them, remove $i from Genomics$i$popname.
		echo -e "Analysis of chr$i complete.\n"
	done
	}

#############################
## ANALYSIS BY POPULATIONS ##
#############################

if [ "$POP" == "FALSE" ]; then 	
	genome_analysis
elif [ "$POP" == "5" ]; then
	# Super-populations (5) are parsed from PopulationIndividualsList.panel expected to be at $gpdat/Others
	superpops=$(cd $gpdat/Others && cut -f3 PopulationIndividualsList.panel | sort -u | grep -v "super_pop" | tr '\n' ' ')
	for spop in $superpops; do
		genome_analysis $spop
	done
elif [ "$POP" == "26" ]; then
	# The 1000 GP populations (26)) are parsed from PopulationIndividualsList.panel expected to be at $gpdat/Others
	allpops=$(cd $gpdat/Others && cut -f2 PopulationIndividualsList.panel | sort -u | grep -v "pop" | tr '\n' ' ')
	for npop in $allpops; do
		genome_analysis $npop
	done	
else
	genome_analysis $POP
fi
