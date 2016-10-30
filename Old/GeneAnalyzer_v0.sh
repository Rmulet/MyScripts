#!/bin/bash

# GENE ANALYZER 0.1 - Evaluates natural selection regimes and population genetics statistics in a gene-restricted manner.

# IMPORTANT: This pipeline assumes the usage of a GFF file with UCSC annotation. Since UCSC data are not provided in this format
# by default, the Perl script UCSC_table2GFF3 is used. Alternatively, it could be downloaded and transformed with gtf2gff3.pl
# The same proces

# Add DOWNLOAD and MASK option

display_usage() { 
	echo -e "\nEvaluates natural selection regimes and population genetics statistics in a gene-restricted manner" 
	echo -e "\nUsage:\n $(basename "$0") I[-msk --mask] [-dl --download]"
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

MASK="pilot"
POP="FALSE"
DL="FALSE"

while [[ $# > 0 ]]
do
	case "$1" in
		-msk|--mask)
		MASK=$(echo "$2" | tr '[:upper:]' '[:lower:]') # $1 has the name, $2 the value
		echo -e "1000 GP mask set to $MASK"
		shift 2
		;;
		-dl|--download)
		DL="$2" # $1 has the name, $2 the value
		echo -e "ALL necessary files will be downloaded: $DL"
		shift 2
		;;
		-pop|--population)
		POP="$2" # $1 has the name, $2 the value
		if [ "$POP" == "26" ]; then echo "All populations will be analysed"
		elif [ "$POP" == "5" ]; then echo "The 5 super-populations will be analysed"
		else echo -e "Population $POP will be analysed"; fi
		shift 2
		;;		
		*) # No more options
	    ;;
	esac
done

##############################
## VARIABLES AND DATA PATHS ##
##############################

gpraw="/home/roger/Documents/2_GenomicsData/1000GP/Chromosomes" # VCF files from 1000 GP divided by chromosomes
gpdat="/home/roger/Documents/2_GenomicsData/1000GP" # No files required 
alnraw="/home/roger/Documents/2_GenomicsData/Alns/Chromosomes" # Human-chimp alignment (MFA.GZ) divided by chromosomes
alndat="/home/roger/Documents/2_GenomicsData/Alns" # Contains FASTA files (FA.GZ/FA)
finaldir="/home/roger/Documents/2_GenomicsData/Final/GeneByGene" # Contains GFF files (can be removed with some tweaking of GFFtoFASTA)

maskdir=$gpdat/Masks/FASTA

downloader() {

wget -nc -nd -r -l 1 -A "ALL.chr*" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

}


#genome_analysis() {
	for i in `seq 1 22` X Y; do
		i=21
		# POLYMORPHISM #
		echo -e "Processing polymorphism data: chr$i" 
		gpfile=chr$i\_gp.vcf.gz # Name of the filtered file
		cd $gpdat
		if [ ! -e chr$i\_gp.vcf.gz ] # Only if the files have not been previously filtered. Not quotes because it's not a variable.
			then
			echo -e "Filtering out inbred individuals in chr$i"
			cd $gpraw
			gpstart=$(ls -d *chr$i.*vcf.gz) # Starting file
			bcftools view -Oz -S ^InbredIndividuals.txt $gpstart > $gpfile # Filter out inbred individuals
			mv $gpfile $gpdat; cd $gpdat
			tabix -p vcf $gpfile
		fi

		if [ "$POP" != "FALSE" ]; then # Note that "PopulationIndividualsList.panel" is assumed to contain individuals and populations
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

		# PREANALYSIS (GFF to FASTA)	
		cd $finaldir
		if [ ! -e "hg19_knownGene.gff3.gz" ]; then # If no annotation is available, then use this script to download and convert to GFF
		ucsc_table2gff3.pl --ftp known --db hg19 --gz
		fi	
		echo -e "Extracting the annotation file of chr$i"
		zcat hg19_knownGene.gff3.gz | grep -P "^chr$i\t" | sed "s/^chr$i/$i/" > chr$i.gff  # Match the ID format of GFF to the VCF files
		cp $alndat/chr$i.fa $finaldir # Copy the FASTA sequence of the chromosome		
		echo -e "Generating the pseudo-FASTA file"
		GFFtoFASTA8.R chr$i.gff chr$i.fa # Convert FASTA to GFF
		
		# MERGE AND ANALYSIS #
		echo -e "Preparing the accessibility mask for chr$i"
		maskfile=$(cd $maskdir && ls -d *chr$i.*)	
		ln -s $gpdat/$gpfile $gpfile;  ln -s $alndat/$alnfile $alnfile # Create symbolic links for the 1000GP and Human-Chimp data
		ln -s $gpdat/$gpfile.tbi $gpfile.tbi;  ln -s $alndat/$alnfile.tbi $alnfile.tbi # Create symbolic links for the index files
		ln -s $maskdir/$maskfile # Create symbolic links for the index files
		echo -e "Analysing polymorphism and divergence in chr$i"
		echo $gpfile
		GeneByGene7.R $gpfile $alnfile $maskfile $i
		echo -e "Analysis of chr$i complete.\n"
# done
	#}

