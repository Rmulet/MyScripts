#!/bin/bash

# VERSION 0.5B - Merges the selected regions of the two VCF files with bcftools.
# WARNING: This script requires a terminating newline. Information on how to circumvent this limitation is available here:
# http://stackoverflow.com/questions/4165135/how-to-use-while-read-bash-to-read-the-last-line-in-a-file-if-there-s-no-new
# WARNING: This script assumes that the alignment file only contains the "Mouse" column, i.e. not the "Human". Otherwise, 
# there will be issues with the merge because only the "Mouse" column is removed in the split.

display_usage() { 
	echo -e "\nThis script merges regions of two compressed VCF files defined in a BED file using BCFtools. To work, it must be provided with one BED file containing positions in the 2nd and 3rd columns, and two VCF files." 
	echo -e "\nUsage:\n $(basename "$0") [BED file] [1000GP VCF.GZ file] [Alignment VCF.GZ file] [-c option]" 
	echo -e "\nOptions:\n -c, -cnvs \t Remove copy number variants (CNV) from the 1000GP file,\n\t\t which may extend beyond the limits of the interval" 
	} 

bedfile=$1
gpfile=$2 # 1000 GP file
alnfile=$3 # Alignment human-mouse file
k=1

# If less than two arguments supplied, display usage 
	if [  $# -le 2 ] 
	then 
		echo -e "\nERROR: Missing arguments"
		display_usage
		exit 1
	fi 

# Check whether user had supplied -h or --help . If yes display usage 
	if [[ ( $# == "--help") ||  $# == "-h" ]] 
	then 
		display_usage
		exit 0
	fi 	

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
	exit
fi

while read chrom pos1 pos2 level
do
	# MERGE THE VCF FILES CORRESPONDING TO THE TWO REGIONS
	echo $k
	echo $chrom $pos1 $pos2
	chrom="${chrom##*[A-Za-z]}" # To extract the chromosome number
	echo $chrom
	if [ "$4" == "-c" ] || [ "$4" == "-cnvs" ]; then
		bcftools merge -r $chrom:$pos1-$pos2 $gpfile $alnfile | grep -v "<CN" > merge.vcf
	else	
		bcftools merge -Ou -o merge.vcf -r $chrom:$pos1-$pos2 $gpfile $alnfile
		echo "1st merge"
	fi
	# REPLACE THE MISSING POSITIONS OF 1000GP WITH 0
	bcftools view -s "^Mouse" merge.vcf > frag.vcf # Remove the Mouse samples
	bcftools +setGT frag.vcf -o filled.vcf.gz -Oz -- -t ./. -n 0  # Replace . with 0
	tabix -f -p vcf filled.vcf.gz
	bcftools merge -r -Oz $chrom:$pos1-$pos2 filled.vcf.gz $alnfile > merge.$k.vcf.gz # Merge filled.vcf with alignment again
	echo "2nd merge"
	((k++))
done < "$bedfile"

# CLEANING AND PREPARING FOR THE NEXT STEP
rm merge.vcf frag.vcf filled.vcf.gz filled.vcf.gz.tbi # Removing temporal files
tabix -p vcf merge.$k.vcf.gz # Indexing for PopGenome

 # bcftools +setGT combined.vcf -o filled.vcf -- -t . -n 0
 # bcftools view -s Mouse merge.1.vcf
 # time tabix -h chr22_aln2.vcf.gz 22:20000000-25000000 | bgzip > test.vcf.gz
 # time cut -f1-2504 test.vcf > a.vcf --> 1'' (but the information of the headers is no longer true)
 # time bcftools view -s "^Mouse" test.vcf.gz > b.vcf --> 5''
