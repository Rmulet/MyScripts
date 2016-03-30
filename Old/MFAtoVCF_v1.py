#!/usr/bin/env python3

# This script converts multi-FASTA alignments (.mfa) files structured as pairs of aligned sequences,
# each of which covers a region of a chromosome. Since the data from the 1000 GP is divided into chromosomes, 
# it is also convenient to have the mouse-human alignments as several VCFs for their subsequent merging.

# WARNING: The tool SNP-sites and all required inputs have to be in the same folder
# WARNING: It is advised to executed this tool on Python3.

import argparse
import subprocess
import re
import os

#######################
## PARSING ARGUMENTS ##
#######################

parser = argparse.ArgumentParser(description='Obtain a VCF file from a MFA using the SNP-sites tool. Currently it is only prepared to handle single .mfa files containing information from one chromosome.')
parser.add_argument('-i','--input',type=str,required=True,help="Input file, in MFA format")
parser.add_argument('-o','--output',type=str,default="output.vcf",help="Output file, in VCF format")
parser.add_argument('-c','--chrom',type=str,default="22",help="Chromosome identifier")

args = parser.parse_args()

############################
## INITIALIZING VARIABLES ##
############################

# VCF MANIPULATION:
trigger = 0
counter = 0
templines = ""
totalrows = ""
# VCF FILE HEADER:
header = ""
version = ""
chrom = args.chrom
seqlength = 0

def vcfprocess(seqlength,totalrows):
	tfile = open("tempin.fasta",'w') # Stores the pair of sequences in a temporal file
	tfile.write(templines) # Writes the temp variable in a .fasta file
	tfile.close()
	subprocess.call(['snp-sites','-v','-o','tempout.vcf','tempin.fasta']) # Calls snp-sites and generates a VCF output
	with open("tempout.vcf") as tempvcf: # Opens the VCF output
		for row in tempvcf:
			length = re.search(r'length=(\d+)>',row)
			if length: 
				seqlength = seqlength + int(length.group(1))
			if row[0]!= "#":
				rsplit = row.split()
				pos = int(header) + int(rsplit[1]) -1 # Arithmetical sum!!! We subtract one because it's a 1-based system
				newrow = '\n' + chrom + '\t' + str(pos) + '\t' + '\t'.join(rsplit[2:]) 
				totalrows = totalrows + newrow					
	return(seqlength,totalrows)

#####################################
## CREATING VCF FILES BY SEQ PAIRS ##
#####################################

with open(args.input,'r') as file: 
	for line in file:
		if line[0]==">":
			trigger+=1
		if trigger==3:
		# TAKES THE TWO FASTA SEQUENCES AND TRANSFORM TO VCF FILE	
			seqlength,totalrows = vcfprocess(seqlength,totalrows)
			trigger=1 # Sets the trigger to 1 again					
			templines = "" # Resets the templines variable to empty
			counter+=1
		# STORES THE POSITION OF THE CURRENT FRAGMENT:
		if trigger == 1: # Stores the length of every aligned region
			match = re.search(r'(chr\d+:)(\d+)-',line)	# I've replaced \d\d with \d+
			if match: header = match.group(2)
		# TEMPORARILY STORES THE CURRENT LINE:
		templines=templines+line
 	# END OF FILE:
	seqlength,totalrows = vcfprocess(seqlength,totalrows)
	counter+=1	
	templines = ""

#################################
## CREATING THE FINAL VCF FILE ##
#################################

out = open(args.output,'w') # Open the final output file. 
out.write('##fileformat=VCFv4.1\n') # We assume it is always going to be the same version (but it can be modified)
out.write('##contig=<ID=%s,length=%d>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' % (args.chrom,seqlength))
out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHuman\tMouse')
out.write(totalrows)
# Deleting temporal files required for SNP-sites
os.remove("tempin.fasta")
os.remove("tempout.vcf")

print("Execution complete. %d FASTA pair(s) processed." %counter)
exit()