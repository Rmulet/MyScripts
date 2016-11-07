#!/bin/bash

# GENE ANALYZER 0.5 - Evaluates natural selection regimes and population genetics statistics in a gene-restricted manner.

# IMPORTANT: This pipeline assumes the usage of a GFF file with UCSC annotation. Since UCSC data are not provided in this format
# by default, the Perl script UCSC_table2GFF3 is used. Alternatively, it could be downloaded and transformed with gtf2gff3.pl
# UPDATE: Symbolic links have been removed. Now the scrips directly access the input from their real location.
# UPDATE: GeneAnalyzer has implemented the CHR and MASK arguments, allowing the user to futher customize the analysis.

# WARNING: Two separate processes of PopGenome CANNOT run in PARALLEL in the same folder.

display_usage() { 
	echo -e "\nEvaluates natural selection regimes and population genetics statistics in a gene-restricted manner" 
	echo -e "\nUsage:\n $(basename "$0") [-chr --chromosome] [-pop --population] [-msk --mask] [-dl --download]"
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
CHR=`seq 1 22` # ALL chromosomes
DB=Gene

while [[ $# -gt 1 ]]
do
	case "$1" in
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
		else echo -e "Population $POP will be analysed"; fi
		shift
		;;
		-chr|--chromosome)
		CHR="$2"
		echo -e "Only chromosome $CHR will be analysed"
		shift
		;;
		*) # No more options
		;;
	esac
shift
done

if [[ -z "$CHR" ]]; then
	echo -e 'Error: You have not specified a chromosome to perform the analysis'
	exit -1
fi

if [[ "$DL" != "TRUE" ]] && [[ "$DL" != "FALSE" ]]; then
	echo -e 'Error: The -dl option must be either TRUE or FALSE'
	exit -1
fi


##############################
## VARIABLES AND DATA PATHS ##
##############################

WORKING="$HOME/Genomics"
 
gpraw="$WORKING/1000GP/Chromosomes" # VCF files from 1000 GP divided by chromosomes
gpdat="$WORKING/1000GP" # No files required 
alnraw="$WORKING/Chromosomes" # Human-chimp alignment (MFA.GZ) divided by chromosomes
alndat="$WORKING/Alns" # Contains FASTA files (FA.GZ/FA)
finaldir="$WORKING/Final/GeneByGene"

# maskdir="$gpdat/Masks/FASTA"

## DOWNLOAD FILES [OPTIONAL]
if [ "$DL" == "TRUE" ]; then
	echo "The files required for the analysis will be downloaded"	
	cd $gpraw
	wget -nc -nd -r -l 1 -A "ALL.chr*" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
	cd $gpdat/Others
	wget -nc ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel -O "PopulationIndividualsList.panel"
	cd $alnraw
	wget -e robots=off -nc -nd -r -l1 -np -A chr??.mfa.gz,chr?.mfa.gz http://pipeline.lbl.gov/data/hg19_panTro4/ # Alignments - VISTA Browser 
	wget -e robots=off -nc -nd -r -l1 -np -A chr??.fa.gz,chr?.fa.gz ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/# Sequence - UCSC
	cd $gpdata/Masks
	wget -nc -nd -r -l0 -np -A strict,pilot ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/
fi

###################
## DATA ANALYSIS ##
###################

cd $finaldir
# touch "Warnings.dat"
START=$(date +%s)

genome_analysis() {
	for i in $1; do # Only autosomal chromosomes
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
			popname=$2
			popins=$(cd $gpdat/Others && grep $popname PopulationIndividualsList.panel | cut -f1 | tr '\n' ',' | sed 's/,$//' )  # List of individuals in that population
			echo "Extracting the individuals of the selected population: $popname"
			if [ ! -e "chr$i$popname.vcf.gz" ]; then
				bcftools view -Oz --force-samples -s $popins chr$i\_gp.vcf.gz > chr$i$popname.vcf.gz # Some have been removed because they are inbred (force-samples to skip)
				tabix -p vcf chr$i$popname.vcf.gz
			fi	
			gpfile=chr$i$popname.vcf.gz # In population mode, then gpfile is the 
			echo $gpfile	
		fi

		 # DIVERGENCE #
		echo -e "Processing divergence data: chr$i" 
		alnfile=chr$i\_aln.vcf.gz # Alignment file in VCF
		echo $alnfile
		cd $alndat
		if [ ! -e "$alnfile" ]; then
			cd $alnraw
			alnstart=chr$i\_aln.mfa # Alignment file in MFA	
			gunzip -c chr$i.mfa.gz | grep -v "score" > $alnstart # Remove 'score' lines from the MFA file
			mv $alnstart $alndat; cd $alndat
			gunzip chr$i.fa.gz # Uncompress the FASTA sequence
			echo -e "Converting .MFA alignment to VCF for chr$i"
			MFAtoVCF.py -s chr$i.fa -q Chimp -c $i $alnstart # Convert MFA to VCF. VCF.GZ file is automatically tabixed.
			rm $alnstart
		fi

		# PREANALYSIS (GFF to FASTA) #	
		cd $finaldir
		echo -e "Extracting the annotation file to a sequence format"
		if [ ! -e "gffseq_chr$i.RData" ]; then
			cp $alndat/chr$i.fa $finaldir # Copy the FASTA sequence of the chromosome
			GFFtoFASTA.R chr$i.fa $i # Puts the GFF annotation in a sequence
		fi

		if [[ "$CHR" != `seq 1 22` ]]; then # To allow external parallelization (multiple instances of Popgenome)
		        mkdir -p $finaldir/chr$i
		        ln -s $finaldir/GenesTable.RData $finaldir/chr$i
		        ln -s $finaldir/gffseq_chr$i.RData $finaldir/chr$i
		        finaldir="$WORKING/Final/GeneByGene/chr$CHR"
			cd $finaldir
		fi

		# MERGE AND ANALYSIS #
		echo -e "Preparing the accessibility mask for chr$i"
		maskfile=$gpdat/Masks/FASTA/$(cd $gpdat/Masks/FASTA && ls -d *chr$i.$MASK*fasta*) # Depends on the chosen criteria. Underscore must be escaped.
		echo -e "Analysing polymorphism and divergence in chr$i"
		echo $gpfile $alnfile ${maskfile##/*} $i GeneData_chr$i $POP  # GPFILE ALNFILE MASK CHROM DB POP
		GeneByGene.R $gpdat/$gpfile $alndat/$alnfile $maskfile $i GeneData_chr$i $popname
		echo -e "Analysis of chr$i complete.\n"		
	done
}

#############################
## ANALYSIS BY POPULATIONS ##
#############################

if [ "$POP" == "FALSE" ]; then  
    genome_analysis $CHR
elif [ "$POP" == "5" ]; then
    # Super-populations (5) are parsed from PopulationIndividualsList.panel expected to be at $gpdat/Others
    superpops=$(cd $gpdat/Others && cut -f3 PopulationIndividualsList.panel | sort -u | grep -v "super_pop" | tr '\n' ' ')
    for spop in $superpops; do
            genome_analysis $CHR $spop
    done
elif [ "$POP" == "26" ]; then
	# The 1000 GP populations (26)) are parsed from PopulationIndividualsList.panel expected to be at $gpdat/Others
	allpops=$(cd $gpdat/Others && cut -f2 PopulationIndividualsList.panel | sort -u | grep -v "pop" | tr '\n' ' ')
	for npop in $allpops; do
			genome_analysis $CHR $npop
	done    
else
	genome_analysis $CHR $POP
fi

END=$(date +%s)
DIFF=$(( END - START ))
echo "Runtime: $DIFF seconds"
