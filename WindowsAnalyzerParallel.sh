#!/bin/bash

# WINDOWS ANALYZER PARALLEL 0.1 - Calls the VCFmerger.sh script for the entire genome in parallel using xargs

# UPDATE 20161025: Feature for specific chromosome analysis added
# UPDATE 20161027: Implemented the MASK feature. Compatibility with Andromeda evaluated and minor issues fixed.

display_usage() {
	echo -e "\nThis script analyses patterns of variation along the entire genome by repeatedly calling VCFmerger.sh for each chromosome"
	echo -e "\nUsage:\n $(basename "$0") [-w --window] [-msk --mask] [-dl --download]"
	echo -e "\nOptions:\n -w, --window \t Specifies the size of the window to be analyzed [10000]"
	echo -e " -msk, --mask \t Indicates what 1000 Genomes Project mask should be used (Pilot/Strict) [Pilot]"
	echo -e " -dl, --download \t Determines whether the necessary files should be downloaded from predefined URLs [FALSE]"
	echo -e " -pop, --population \t Conducts the analysis by population (26/5/name/FALSE) [FALSE]"
	echo -e " -chr, --chromosome \t Analyses a single chromosome (in format NN, e.g. 22) [ALL]"
	echo -e " -t, --threads \t Number of chromosomes analysed in parallel [2]"
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

export WINDOW=10000
export MASK="pilot"
export POP="FALSE"
export DL="FALSE"
export CHR="`seq 1 22` X Y" # ALL chromosomes
export THREADS=2

while [[ $# -gt 1 ]] # next two arguments (window + size)

do
	case "$1" in
		-w|--window)
		export WINDOW="$2" # $1 has the name, $2 the value
		echo -e "Window size set to $WINDOW"
		shift # next two arguments (window + size)
		;;
		-msk|--mask)
		export MASK=$(echo "$2" | tr '[:upper:]' '[:lower:]') # $1 has the name, $2 the value
		echo -e "1000 GP mask set to $MASK"
		shift
		;;
		-dl|--download)
		DL="$2" # $1 has the name, $2 the value
		echo -e "ALL necessary files will be downloaded: $DL"
		shift
		;;
		-pop|--population)
		export POP="$2" # $1 has the name, $2 the value
		if [ "$POP" == "26" ]; then echo -e "All populations will be analysed\n"
		elif [ "$POP" == "5" ]; then echo -e "The 5 super-populations will be analysed\n"
		else echo -e "Population $POP will be analysed\n"; fi
		shift
		;;		
		-chr|--chromosome)
		export CHR="$2" # $1 has the name, $2 the value
		echo -e "Only chromosome $CHR will be analysed\n"
		shift
		;;
		-t|--threads)
		THREADS="$2"
		echo -e "The analysis will be performed on $THREADS chromosomes at once"
		;;
		*) # No more options
		;;
	esac
shift
done

sleep 5 # Delay for 5 seconds

# Redirect stdout ( > ) into a named pipe ( >() ) running "tee".
rm -f WindowsAnalyzer.log
exec > >(tee -a WindowsAnalyzer.log) 2>&1

##############################
## VARIABLES AND DATA PATHS ##
##############################

export WORKING="$HOME/Genomics"

export gpraw="$WORKING/1000GP/Chromosomes" # VCF files from 1000 GP divided by chromosomes
export gpdat="$WORKING/1000GP" # No files required 
export alnraw="$WORKING/Alns/Chromosomes" # Human-chimp alignment (MFA.GZ) divided by chromosomes
export alndat="$WORKING/Alns" # Contains FASTA files (FA.GZ/FA)
export finaldir="$WORKING/Final" # Contains GFF files

export BCFTOOLS="/home/roger/Software/bcftools"

export maskfile=$gpdat/Masks/$(cd $gpdat/Masks/ && ls -d *$MASK\_mask.whole_genome.bed) # Depends on the chosen criteria. 

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
	i=$1
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
		$BCFTOOLS/bcftools view -Oz --force-samples -s $fem $gpfile > chr$i.temp.vcf.gz # Remove female individuals
		mv chr$i.temp.vcf.gz $gpfile
		tabix -p vcf chr$i.temp.vcf f $gpfile
	elif [ "$i" == "Y" ]; then
		echo -e "Converting chromosome Y to pseudo-diploid"
		zcat $gpfile | perl -p -e 's/\t([01.])(?=[\n|\t])/\t\1\|\1/g' > chr$i.temp.vcf # Duplicate haploid individuals
		bgzip chr$i.temp.vcf; mv chr$i.temp.vcf.gz $gpfile
		tabix -p vcf $gpfile
	fi

	if [ "$POP" != "FALSE" ]; then # Note that "PopulationIndividualsList.panel" is assumed to integrated_call_samples_v3ontain
		popname=$2
		echo "1"
		popins=$(cd $gpdat/Others && grep $popname PopulationIndividualsList.panel | cut -f1 | tr '\n' ',' | sed 's/,$//' )  # List of individuals in that population
		echo "Extracting the individuals of the selected population: $popname"
		$BCFTOOLS/bcftools view -Oz --force-samples -s $popins chr$i\_gp.vcf.gz > chr$i$popname.vcf.gz # Some have been removed because they are inbred (force-samples to skip)
		gpfile=chr$i$popname.vcf.gz # In population mode, then gpfile is the 
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
	if [ ! -e "gffseq_chr$i.RData" ]; then
		cp $alndat/chr$i.fa $finaldir # Copy the FASTA sequence of the chromosome		
		GFFtoFASTA.R chr$i.fa $i # Puts the GFF annotation in a sequence
	fi			

	# MERGE AND ANALYSIS #
	echo -e "Preparing the accessibility mask for chr$i"
	grep "chr$i" $maskfile > $MASK\_mask.chr$i.bed # Extract the chromosome of interest from the mask.		
	echo -e "Analysing polymorphism and divergence in chr$i"

	mkdir -p $finaldir/chr$i
        ln -s $finaldir/gffseq_chr$i.RData $finaldir/chr$i # Link to the GFF file
	ln -s $finaldir/$MASK\_mask.chr$i.bed $finaldir/chr$i # Link to the mask file
        finaldir="$WORKING/Final/chr$CHR"
	
	cd $finaldir               
	VCFmerger.sh $MASK\_mask.chr$i.bed $gpdat/$gpfile $alndat/$alnfile -w $WINDOW -db Genomics$i$popname # Merges and analyzes variation data

	if [[ $? -ne 0 ]]; then # Stop the execution of the script if VCFmerger fails
		echo 'Error: VCFmerger.sh failed!'
		exit -1
	fi

	# The $popname variable will be ommitted if it is not declared (i.e. analysis conducted for 2,504 individuals)
	# By default, the results will be stored in separate tables. To merge them, remove $i from Genomics$i$popname.
	echo -e "Analysis of chr$i complete.\n"
}

export -f genome_analysis

#############################
## ANALYSIS BY POPULATIONS ##
#############################

if [ "$POP" == "FALSE" ]; then
	echo $CHR | xargs -n 1 -P$THREADS bash -c 'genome_analysis "$@"' _
elif [ "$POP" == "5" ]; then
	# Super-populations (5) are parsed from PopulationIndividualsList.panel expected to be at $gpdat/Others
	superpops=$(cd $gpdat/Others && cut -f3 PopulationIndividualsList.panel | sort -u | grep -v "super_pop" | tr '\n' ' ')
	for spop in $superpops; do
		echo $CHR | xargs -n 1 -P$THREADS bash -c 'genome_analysis "$@" $spop' _
	done
elif [ "$POP" == "26" ]; then
	# The 1000 GP populations (26)) are parsed from PopulationIndividualsList.panel expected to be at $gpdat/Others
	allpops=$(cd $gpdat/Others && cut -f2 PopulationIndividualsList.panel | sort -u | grep -v "pop" | tr '\n' ' ')
	for npop in $allpops; do
		echo $CHR | xargs -n 1 -P$THREADS bash -c 'genome_analysis "$@" $npop' _
	done
else
	echo $CHR | xargs -n 1 -P$THREADS bash -c 'genome_analysis "$@" $POP' _
fi
