#!/bin/bash

#--------------------------------------------------------------------------------------------------
#Get input arguments
#--------------------------------------------------------------------------------------------------

display_usage() {

cat << HEREDOC

	This script analyses patterns of variation along the entire genome by repeatedly calling VCFmerger.sh for each chromosome

	Usage: exact_windows.sh [-w --window] [-msk --mask] [-dl --download]

	Options: 
	-w, --window			Specifies the size of the window to be analyzed [10000]
	-msk, --mask			Indicates what 1000 Genomes Project mask should be used (Pilot/Strict) [Pilot]
	-dl, --download 		Determines whether the necessary files should be downloaded from predefined URLs [FALSE]
	-pop, --population		Conducts the analysis by population (26/5/name/ALL) [ALL]
	-chr, --chromosome		Analyses a single chromosome (in format NN, e.g. 22) [ALL]
	-db, --database			Specifies the name of the data table where the results will be stored [GenomicsWhole]"

HEREDOC

}

WINDOW=10000
MASK="pilot"
POP="ALL"
DL="FALSE"
CHR="`seq 1 22` X Y" # ALL chromosomes
MKT=TRUE # MKT will be calculated by default
DB="GenomicsWhole"

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
		-h|--help)
		display_usage
		exit
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

######################
## SCRIPT EXECUTION ##
######################

## Generate a fake VCF file

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

		# GENERATION OF MASKED POSITIONS FILE

		MASK_POS=$gpdat/Masks/POS/chr$i.pos_notac.vcf.gz
		if [[ ! -e $MASK_POS ]]; then
		
			NSAMPLES=$(bcftools query -l $gpdat/$gpfile | wc -l)	
			SIZE=$(grep chr$i $gpdat/Masks/hg19.chrom.sizes | cut -f2)
			START=$(grep chr$i $maskfile | head -1 | cut -f2)

			# awk -v var="chr$i" '$1==var' $maskfile > $MASK\_mask.chr$i.bed # Extract the chromosome of interest from the mask
			echo -e "Extracting masked coordinates in chr$i"
			echo -e "chr$i\t$START\t$SIZE" > chr${i}_total.bed
			bedtools subtract -a chr${i}_total.bed -b $maskfile > chr${i}_notac.bed

			time awk -v var="chr$i" '$1==var{S=int($2)+1;E=($3);for(i=S;i<=E;++i) printf("%s %d %d\n",$1,i,i);}' chr${i}_notac.bed | sort | uniq | sort -t ' ' -k1,1 -k2,2n > $gpdat/Masks/POS/chr$i.pos_notac.bed
			rm chr${i}_total.bed chr${i}_notac.bed

			echo -e "Generating VCF file with masked positions of chr$i"
			echo "##fileformat=VCFv4.1\n" > $gpdat/Masks/POS/chr$i.pos_notac.vcf
			echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"' >> $gpdat/Masks/POS/chr$i.pos_notac.vcf
			bcftools view -h $gpdat/$gpfile | grep -v '##' >> $gpdat/Masks/POS/chr$i.pos_notac.vcf # Create the header of the file
			time cat $gpdat/Masks/POS/chr$i.pos_notac.bed | while read -a P; do echo -e -n "${P[0]}\t${P[1]}\t.\t" && samtools faidx /ucsc/hg19/hg19.fa "${P[0]}:${P[1]}-${P[1]}"\
			| grep -v '>' | awk -v N=$NSAMPLES 'BEGIN{NSAMPLES=N;} {printf("%s\tN\t.\t.\t.\tGT",$0); for(i=1;i<=NSAMPLES;i++) printf("\t./."); printf("\n");}'; done >> $gpdat/Masks/POS/chr$i.pos_notac.vcf

			bgzip $gpdat/Masks/POS/chr${i}.pos_notac.vcf
			tabix -p vcf $MASK_POS
			
		fi

		# FILTERING SEXUAL CHROMOSOMES

		if [ "$i" == "X" ]; then # Remove MALES from the X chromosome file
			echo -e "Excluding males from chromosome X"
			fem=$(cd $gpdat/Others && grep "female" PopulationIndividualsList.panel | cut -f1 | tr '\n' ',')
			echo $fem
			$BCFTOOLS/bcftools view -Oz --force-samples -s $fem $gpfile > chr$i.temp.vcf.gz # Remove female individuals
			mv chr$i.temp.vcf.gz $gpfile
			tabix -p vcf $gpfile

			$BCFTOOLS/bcftools view -Oz --force-samples -s $fem $MASK_POS > $gpdat/Masks/POS/chr$i.pos_temp.vcf.gz # Remove female individuals
			mv $gpdat/Masks/POS/chr$i.pos_temp.vcf.gz $MASK_POS
			tabix -p vcf $MASK_POS

		elif [[ "$i" == "Y" ]] && [[ $(zcat $gpfile | grep -v '#' | head -1 | awk '{if ($10 ~ /[0-4.]\|[0-4.]/) print "DIPLOID"; else print "HAPLOID"}') == "HAPLOID" ]]; then
			echo -e "Converting chromosome Y to pseudo-diploid"
			zcat $gpfile | awk '$5 !~ "<CN"' |  perl -pe 's/\t([0-4.])(?=[\n|\t])/\t\1\|\1/g' > chr$i.temp.vcf # Duplicate haploid individuals
			mv $gpfile chr$i.haploid_gp.vcf.gz
			bgzip chr$i.temp.vcf; mv chr$i.temp.vcf.gz $gpfile
			tabix -p vcf $gpfile
		fi

		# FILTERING BY POPULATION

		if [ "$POP" != "ALL" ]; then # Note that "PopulationIndividualsList.panel" is assumed to integrated_call_samples_v3ontain
			popname=$1
			popins=$(cd $gpdat/Others && grep $popname PopulationIndividualsList.panel | cut -f1 | tr '\n' ',' | sed 's/,$//' )  # List of individuals in that population

			if [ ! -e "chr$i$popname.vcf.gz" ]; then		
				echo "Extracting the individuals from $gpfile of the selected population: $popname"; sleep 2
				$BCFTOOLS/bcftools view -Oz --force-samples -s $popins chr$i\_gp.vcf.gz > chr$i$popname.vcf.gz # Some have been removed because they are inbred (force-samples to skip)
				tabix -p vcf chr$i$popname.vcf.gz

			fi

			if [ ! -e $gpdat/Masks/POS/chr$i${popname}_notac.vcf.gz ]; then
				echo "Extracting the individuals from $MASK_FILE of the selected population: $popname"; sleep 2
				$BCFTOOLS/bcftools view -Oz --force-samples -s $popins $MASK_POS > $gpdat/Masks/POS/chr$i${popname}_notac.vcf.gz
				tabix -p vcf $gpdat/Masks/POS/chr$i${popname}_notac.vcf.gz

			fi

			gpfile=chr$i$popname.vcf.gz # In population mode, then gpfile is the file generated for that population
			MASK_FILE=$gpdat/Masks/POS/chr$i${popname}_notac.vcf.gz
			echo $gpfile	

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

        if [[ "$CHR" != `seq 1 22` ]]; then # To allow external parallelization IF only one chr is analyzed
            mkdir -p $finaldir/chr$i
            ln -s $finaldir/gffseq_chr$i.RData $finaldir/chr$i # Link to the GFF file
            finaldir="$WORKING/Final/chr$CHR"
			cd $finaldir
        fi
		
#		bedtools intersect -v -a $gpfile -b $gpdat/Masks/POS/chr$i.pos_masked.bed > ${gpfile%%.vcf.gz}_masked.vcf.gz # Remove masked positions
#		gpfile=${gpfile%%.vcf.gz}_masked.vcf.gz

		SIZE=$(grep chr$i $gpdat/Masks/hg19.chrom.sizes | cut -f2)

		for (( pos1=1 ; pos1 <= $SIZE ; pos1=pos1+$WINDOW )); do

			# MERGE THE VCF FILES, REPLACING MISSING GENOTYPES
			# The genome accessibility mask is in exclusive 1-based format, i.e. the last position is not included [half open]. 
			# Since bcftools assumes inclusive 1-based format [closed], we subtract 1 from the end coordinate.
			pos2=$((pos1+$WINDOW-1))
			chrom=$i
			bcftools merge -Oz --missing-to-ref -o merge.$pos1.vcf.gz -r $chrom:$pos1-$pos2 $gpdat/$gpfile $alndat/$alnfile

			if [[ $? -ne 0 ]]; then
				echo 'Error: Bcftools merge failed!'
				exit -1
			else
				echo -e "Files merged: merge.$pos1.vcf.gz generated"
			fi
	
			tabix -p vcf merge.$pos1.vcf.gz # Tabixing for analysis with PopGenome

			# R ANALYSIS OF NUCLEOTIDE VARIATION:
			# PopGenome does not use the first position [left open], so we subtract -1 from its initial position (internal).

			echo merge.$pos1.vcf.gz $chrom $pos1 $pos2 $WINDOW $DB $MKT $POP
			VCFAnalysis.R merge.$pos1.vcf.gz $chrom $pos1 $pos2 $WINDOW $DB $MKT $POP  # Avoid the message visualization!! 

		done

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


