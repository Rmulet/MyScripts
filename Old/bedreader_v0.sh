#!/bin/bash

# WARNING: This script requires a terminating newkline. Information on how to circumvent this limitation is available here:
# http://stackoverflow.com/questions/4165135/how-to-use-while-read-bash-to-read-the-last-line-in-a-file-if-there-s-no-new

display_usage() { 
	echo -e "\nThis script splits VCF files into regions of interest and creates tabixed VCF files. To work, it must be provided with one BED file containing positions in the 2nd and 3rd columns, and two VCF files." 
	echo -e "\nUsage:\n$0 [BED file] [1000 GP VCF file] [Alignment VCF file] [-c option] \n" 
	} 

bedfile=$1
gpfile=$(basename $2 .vcf.gz)
echo $gpfile
alnfile=$(basename $3 .vcf.gz)
echo $alnfile
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
	if [ "$4" == "-c" ] || [ "$4" == "-cnvs" ]; then
		tabix $2 $chrom:$pos1-$pos2 | grep -v "<CN" > $gpfile.$k.vcf # We remove CNVs
		tabix $3 $chrom:$pos1-$pos2 > $alnfile.$k.vcf # We remove CNVs
	else	
		tabix $2 $chrom:$pos1-$pos2 > $gpfile.$k.vcf	
		tabix $3 $chrom:$pos1-$pos2 > $alnfile.$k.vcf 
	fi
	bgzip $gpfile.$k.vcf; bgzip $alnfile.$k.vcf
	tabix -p vcf $gpfile.$k.vcf.gz; tabix -p vcf $alnfile.$k.vcf.gz
	((k++))
	# echo $chrom $pos1 $pos2
done < "$bedfile"


