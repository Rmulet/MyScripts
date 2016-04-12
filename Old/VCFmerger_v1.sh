#!/bin/bash

# VERSION 0.4A - Divides the VCF files into pieces and merges them together.
# WARNING: This script requires a terminating newkline. Information on how to circumvent this limitation is available here:
# http://stackoverflow.com/questions/4165135/how-to-use-while-read-bash-to-read-the-last-line-in-a-file-if-there-s-no-new

display_usage() { 
	echo -e "\nThis script splits VCF files into regions of interest and creates tabixed VCF files. To work, it must be provided with one BED file containing positions in the 2nd and 3rd columns, and two VCF files." 
	echo -e "\nUsage:\n$0 [BED file] [1000 GP VCF file] [Alignment VCF file] [-c option] \n" 
	} 

bedfile=$1
gpfile=$(basename $2 .vcf.gz)
alnfile=$(basename $3 .vcf.gz)
k=1

# if less than two arguments supplied, display usage 
	if [  $# -le 2 ] 
	then 
		display_usage
		exit 1
	fi 

# check whether user had supplied -h or --help . If yes display usage 
	if [[ ( $# == "--help") ||  $# == "-h" ]] 
	then 
		display_usage
		exit 0
	fi 	

while read chrom pos1 pos2 level
do
	echo $k
	echo $pos1 $pos2
	chrom="${chrom##*[A-Za-z]}" # To extract the chromosome number
	echo $chrom
	if [ "$4" == "-c" ] || [ "$4" == "-cnvs" ]; then
		tabix -h $2 $chrom:$pos1-$pos2 | grep -v "<CN" > $gpfile.temp.vcf # We remove CNVs
		tabix -h $3 $chrom:$pos1-$pos2 > $alnfile.temp.vcf # We remove CNVs
	else	
		tabix -h $2 $chrom:$pos1-$pos2 > $gpfile.temp.vcf	
		tabix -h $3 $chrom:$pos1-$pos2 > $alnfile.temp.vcf 
	fi
	bgzip $gpfile.temp.vcf; bgzip $alnfile.temp.vcf
	tabix -p vcf $gpfile.temp.vcf.gz; tabix -p vcf $alnfile.temp.vcf.gz
	vcf-merge $gpfile.temp.vcf.gz $alnfile.temp.vcf.gz > merge.$k.vcf
	rm $gpfile.temp.vcf.gz $alnfile.temp.vcf.gz $gpfile.temp.vcf.gz.tbi $alnfile.temp.vcf.gz.tbi
	echo $k
	((k++))
	# echo $chrom $pos1 $pos2
done < "$bedfile"


