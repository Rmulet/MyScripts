#!/bin/bash

display_usage() {
	echo -e "\nThis script analyses patterns of variation along the entire genome by repeatedly calling VCFmerger.sh for each chromosome"
	echo -e "\nUsage:\n $(basename "$0") [-w --window] [-msk --mask] [-dl --download]"
	echo -e "\nOptions:\n -w, --window \t Specifies the size of the window to be analyzed [10000]"
	echo -e " -msk, --mask \t Indicates what 1000 Genomes Project mask should be used (Pilot/Strict) [Pilot]"
	echo -e " -dl, --download \t Determines whether the necessary files should be downloaded from predefined URLs [FALSE]"
	echo -e " -pop, --population \t Conducts the analysis by population (26/5/name/FALSE) [FALSE]"
	echo -e " -chr, --chromosome \t Analyses a single chromosome (in format NN, e.g. 22) [ALL]"
	}

WINDOW=10000
MASK="pilot"
POP="FALSE"
DL="FALSE"
CHR="`seq 1 22` X Y" # ALL chromosomes

while [[ $# -gt 1 ]] # next two arguments (window + size)

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

WORKING="$HOME/Genomics"

gpraw="$WORKING/1000GP/Chromosomes" # VCF files from 1000 GP divided by chromosomes
gpdat="$WORKING/1000GP" # No files required 
alnraw="$WORKING/Alns/Chromosomes" # Human-chimp alignment (MFA.GZ) divided by chromosomes
alndat="$WORKING/Alns" # Contains FASTA files (FA.GZ/FA)
finaldir="$WORKING/Final" # Contains GFF files

BCFTOOLS="/home/roger/Software/bcftools"

maskfile=$gpdat/Masks/$(cd $gpdat/Masks/ && ls -d *$MASK\_mask.whole_genome.bed) # Depends on the chosen criteria. Underscore must be escaped.
echo $maskfile

######################
## SCRIPT EXECUTION ##
######################

## Generate a fake VCF file

genome_analysis() {
	for i in $CHR; do

		if [[ ! -e $gpdat/Masks/POS/chr$i.pos_masked.vcf ]]; then

			#awk -v var="chr$i" '$1==var' $maskfile > $MASK\_mask.chr$i.bed # Extract the chromosome of interest from the mask
			echo -e "Extracting masked coordinates in chr$i"
			awk -v var="chr$i" '$1==var{S=int($2)+1;E=($3);for(i=S;i<=E;++i) printf("%s %d %d\n",$1,i,i);}' $maskfile | sort | uniq | sort -t ' ' -k1,1 -k2,2n > $gpdat/Masks/POS/chr$i.pos_masked.bed

			echo -e "Generating VCF file with masked positions of chr$i"
			cat $gpdat/Masks/POS/chr$i.pos_masked.bed | while read -a P; do echo -e -n "${P[0]}\t${P[1]}\t.\t" && samtools faidx /ucsc/hg19/hg19.fa "${P[0]}:${P[1]}-${P[1]}" | grep -v '>' | awk 'BEGIN{NSAMPLES=2;} {printf("%s\tN\t.\t.\t.\tGT",$0); for(i=1;i<=NSAMPLES;i++) printf("\t./."); printf("\n");}'; done >> $gpdat/Masks/POS/chr$i.pos_masked.vcf
			
		fi

	done

}

genome_analysis


