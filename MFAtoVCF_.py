#!/usr/bin/env python3

# MFAtoVCF.py, Version 1.5

# This script converts multi-FASTA alignments (.mfa) files structured as pairs of aligned sequences,
# each of which covers a region of a chromosome. Since the data from the 1000 GP is divided into chromosomes, 
# it is also convenient to have the mouse-human alignments as several VCFs for their subsequent merging.

# WARNING: The tool SNP-sites and all required inputs have to be in the same folder
# WARNING: It is advised to executed this tool on Python3.
# WARNING: This version detects overlaps and changes repeats to Rs.

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
parser.add_argument('-t','--test',action='store_true',help="Testing mode, keeps temporal files for examination")

args = parser.parse_args()

############################
## INITIALIZING VARIABLES ##
############################

# VCF MANIPULATION:
trigger = 0 # Trigger=0 - Start; Trigger=1 - 1st human seq; Trigger=2 - 1st mouse seq; Trigger=3 - 2nd human seq
counter = 0
human1,mouse1 = "",""
totalrows = ""
end0 = 0 # Longest end so far
aligned = []
# VCF FILE HEADER:
startpos = ""
headers = [] # First human, second mouse
version = ""
chrom = args.chrom
seqlength = 0

##########################
## FUNCTION DEFINITION ##
#########################

def fastainput(human1,mouse1,end0):
	human2 = mouse2 = mouse3 = ""
	for i,base in enumerate(human1): # REMOVE GAPS IN HUMAN - LENGTH EQUAL TO REGION WIDTH
		if base != "-":
			mouse2 = mouse2 + mouse1[i]
			human2 = human2 + base
	if start < end0: # OVERLAP --> It only works for the NEXT region!
		for i,base in enumerate(mouse2): 
			if base == "-":	
				mouse3 = mouse3 + "X" # X for Gap -- Identify gaps in mice to be removed later
			elif i < (end0-start+1) and base == human2[i]: # Only othe overlapping part of the region
				mouse3 = mouse3 + "R" # P for Polymorphism -- a) Multiple Ps: retain human; b) One P: remove site
			else:
				mouse3 = mouse3 + base	
	else: # NOT OVERLAP
		for i,base in enumerate(mouse2):
			if base == "-":	
				mouse3 = mouse3 + "X" # X for Gap
			else:
				mouse3 = mouse3 + base	
	aligned.append([start,end1])
	if end1 > end0:
		end0 = end1			
	tfile = open("tempin.fasta",'w') # Stores the pair of sequences in a temporal file
	tfile.write('%s%s%s%s'%(headers[0],human2,headers[1],mouse3)) # Writes the temp variable in a .fasta file
	tfile.close()
	if args.test:
		print (headers[0].strip()+':',len(human2))
	print ("Human"+human2+'\n'+"Mouse"+mouse3)
	return(end0,aligned)	

def vcfprocess(seqlength,totalrows):
	subprocess.call(['snp-sites','-v','-o','tempout.vcf','tempin.fasta']) # Calls snp-sites and generates a VCF output
	with open("tempout.vcf") as tempvcf: # Opens the VCF output
		for row in tempvcf:
			length = re.search(r'length=(\d+)>',row)
			if length: 
				seqlength = seqlength + int(length.group(1))
			if row[0]!= "#":
				rsplit = row.split()
				pos = start + int(rsplit[1]) -1 # Arithmetical sum!!! We subtract one because it's a 1-based system
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
		if trigger==3: # Every Human (or reference species) line after the first
		# REMOVE THE GAPS FROM THE HUMAN SEQUENCE:
			end0,aligned = fastainput(human1,mouse1,end0) # Generates the .fasta file 
			seqlength,totalrows = vcfprocess(seqlength,totalrows)
			trigger=1 # Sets the trigger to 1 again					
			human1 = mouse1 = ""
			headers = []
			counter+=1
		# STORES THE POSITION OF THE CURRENT FRAGMENT:
		if trigger == 1 and line[0]==">": # Stores the length of every aligned region
			match = re.search(r'(chr\d+:)(\d+)-(\d+)',line)	# I've replaced \d\d with \d+
			if match: 
				start = int(match.group(2))
				end1 = int(match.group(3))
		# STORE INFORMATION TO CONSTRUCT THE FASTA FILE:	
			headers.append(line)
		elif trigger == 2 and line[0]==">":
			headers.append('\n'+line)
		elif trigger == 1 and line[0] != ">":
			human1=human1+line.strip()
		elif trigger == 2 and line[0] != ">":
			mouse1=mouse1+line.strip()
 	# UPON REACHING THE END OF THE FULE:
	fastainput(human1,mouse1,end0)
	seqlength,totalrows = vcfprocess(seqlength,totalrows)		
	human1 = mouse1 = ""
	headers = []
	counter+=1	
	
#################################
## CREATING THE FINAL VCF FILE ##
#################################

out = open(args.output,'w') # Open the final output file. 
out.write('##fileformat=VCFv4.1\n') # We assume it is always going to be the same version (but it can be modified)
out.write('##contig=<ID=%s,length=%d>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' % (args.chrom,seqlength))
out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHuman\tMouse')
out.write(totalrows)
out.close()
# Deleting temporal files required for SNP-sites
if not args.test:
	os.remove("tempin.fasta")
	os.remove("tempout.vcf")

print("Execution complete. %d FASTA pair(s) processed covering a total of %s nucleotides." % (counter,seqlength))
exit()