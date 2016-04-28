#!/bin/bash

# VERSION 0.6 - Merges the selected regions of the two VCF files with bcftools. Calculates statistics with R.
# WARNING: This script requires a terminating newline in the BED file. Information on how to circumvent this limitation is available here:
# http://stackoverflow.com/questions/4165135/how-to-use-while-read-bash-to-read-the-last-line-in-a-file-if-there-s-no-new
# WARNING: This script requires the usage of the latest experimental version of BCFTOOLS, available at http://pd3.github.io/bcftools/
# WARNING: RegionAnalysis.R is assumed to be in an environment variable. If not, modify the address.

# UPDATE: Argument parsing has been improved. Window size = 200 by default.

# Makes sure that that no repeated rows are added
mysql --user="root" --password="RM333" --database="PEGH" --execute="DROP table IF EXISTS Genomics";

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
	if [ "$(($pos2 - $pos1))" -le "$WINDOW" ]; then # If size region < size window, there is no need to examine it 
		echo -e "$pos2 - $pos1 is inferior to $WINDOW. This segment will be skipped."
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

	# R PROCESSING:
	RegionAnalysis.R merge.$k.vcf.gz $pos1 $pos2 $WINDOW >/dev/null # Avoid the message visualizatio
	# Filename = merge.$k.vcf; ini = pos1; end = pos2 // --slave >/dev/null  [slave to cut startup messages]
	echo -e "Analysis with R complete. Data exported to the MySQL database."
	rm merge.$k.vcf.gz merge.$k.vcf.gz.tbi # Removes the current merge file in order to free disk space

	((k++))
done < "$bedfile"

echo -e "\nAll regions in the BED file merged and analyzed."