#!/bin/bash

# VERSION 0.5B - Merges the selected regions of the two VCF files with bcftools.
# WARNING: This script requires a terminating newline. Information on how to circumvent this limitation is available here:
# http://stackoverflow.com/questions/4165135/how-to-use-while-read-bash-to-read-the-last-line-in-a-file-if-there-s-no-new
# WARNING: This script assumes that the alignment file only contains the "Mouse" column, i.e. not the "Human". Otherwise, 
# there will be issues with the merge because only the "Mouse" column is removed in the split.
# WARNING: This script requires the usage of the latest experimental version of BCFTOOLS, available at http://pd3.github.io/bcftools/
# UPDATE: Argument parsing has been improved.

display_usage() { 
	echo -e "\nThis script merges regions of two compressed VCF files defined in a BED file using BCFtools. To work, it must be provided with one BED file containing positions in the 2nd and 3rd columns, and two VCF files." 
	echo -e "\nUsage:\n $(basename "$0") [BED file] [1000GP VCF.GZ file] [Alignment VCF.GZ file] [-w --window] [-c --cnvs]" 
	echo -e "\nOptions:\n -c, --cnvs \t Remove copy number variants (CNV) from the 1000GP file,\n\t\t which may extend beyond the limits of the interval [False]" 
	echo -e " -w, --window \t Specifies the size of the window to be analyzed [200]"
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

# If less than two arguments supplied, display usage 
if [  $# -le 2 ] 
then 
	echo -e "\nERROR: Missing arguments"
	display_usage
	exit 1
fi 

## ARGUMENT PARSING ##

# POSITIONAL ARGUMENTS:

bedfile=$1; shift # BED file containing the regions
gpfile=$1; shift # 1000 GP file
alnfile=$1; shift # Alignment human-mouse file

# Check whether the supplied files exist:

if [ ! -e "$bedfile" ] || [ ! -e "$gpfile" ] || [ ! -e "$alnfile" ] 
	then
	if [ ! -e "$bedfile" ]
		then
		echo -e "ERROR: $bedfile not found.\n"
	fi
	if [ ! -e "$gpfile" ] 
		then
		echo -e "ERROR: $gpfile not found.\n"
	fi
	if [ ! -e "$alnfile" ] 
	then
	echo -e "ERROR: $alnfile not found.\n"		
	fi
	exit 1
fi

# OPTIONAL ARGUMENTS:

WINDOW=200

while [[ $# > 0 ]]
do
	case "$1" in
		-w|--window)
		WINDOW="$2" # $1 has the name, $2 the value
		echo "Window size set to $WINDOW"
		shift 2 # next two arguments (window + size)
		;;
		-c|--cnvs)
		CNVS="CNVS" # $1 has the name, $2 the value
		echo -e "CNVS will be removed"
		shift # next argument
		;;
		*) # No more options
	    ;;
	esac
done

###############################
## VCF MERGER AND R ANALYSIS ##
###############################

k=1
while read chrom pos1 pos2 level
do
	# EXTRACT BED INFORMATION
	chrom="${chrom##*[A-Za-z]}" # To extract the chromosome number
	echo $k
	echo $chrom $pos1 $pos2 
	if [ "$(($pos2 - $pos1))" -le "200" ]; then # If size region < size window, there is no need to examine it 
		expr $pos2 - $pos1
		continue
	fi	
	# MERGE THE VCF FILES, REPLACING MISSING GENOTYPES
	if [ "$CNVS" == "CNVS" ] ; then
		bcftools merge -Ov --missing-to-ref -r $chrom:$pos1-$pos2 $gpfile $alnfile | grep -v "<CN" > merge.$k.vcf
		bgzip merge.$k.vcf
	else	
		bcftools merge -Oz --missing-to-ref -o merge.$k.vcf.gz -r $chrom:$pos1-$pos2 $gpfile $alnfile
	fi
	echo -e "Files merged: merge.$k.vcf.gz generated"
	tabix -p vcf merge.$k.vcf.gz # Tabixing for analysis with PopGenome
	((k++))

	# R PROCESSING:

done < "$bedfile"

 # TIMINGS: 1 MB region #

 # Merge BCF vs VCF: 1'10'' vs 2'26'' (the difference is less if we output as VCF)
 # Remove mouse BCF vs VCF: 5'36'' vs 8'5''
 # Replace BCF vs VCF: 4'40'' vs 4'30''
 # Index BCF vs VCF: 21'' both
 # Final merge BCF vs VCF: 3'3'' vs 4'24''

 # Total pipeline VCF (Andromeda): 12'19''
 # Total pipeline BCF (Andromeda): 10'28''

 # Script Perl (from VCF, Andromeda): 1'26'' (merge) + 1'18'' (correct) => 2'44'' + 1'6'' (bgzip) + 23'' (tabix) => 4'13''
 # Script Perl (from BCF, Andromeda): 1'18'' (merge) + 1'18'' (correct) => 2'36'' + 1'6'' (bgzip) + 23'' (tabix) => 4'6''

 # Bcftools merge --missing-to-ref (from VCF, not Andromeda): 2'31''
 # Bcftools merge --missing-to-ref (from BCF, Andromeda): 1'44''

 # Bcftools merge --missing-to-ref (from VCF,Andromeda) + bgzip + tabix => 2'56''
 # Bcftools merge --missing-to-ref (from VCF,Andromeda) + bgzip + tabix => 3'8''