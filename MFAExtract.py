#!/usr/bin/env python3

# MFAtoVCF.py, Version 0.0

# This script converts multi-FASTA alignments (.mfa) files structured as pairs of aligned sequences,
# each of which covers a region of a chromosome. Since the data from the 1000 GP is divided into chromosomes, 
# it is also convenient to have the mouse-human alignments as several VCFs for their subsequent merging.

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
human1,mouse1 = "",""
totalrows = ""
human0 = ""
mouse0 = ""

# VCF FILE HEADER:
headers = [] # First human, second mouse
seqlength = ""
chrom = args.chrom
counter = "" # Count the number of processed FASTA regions

#############################
## SNP EXTRACTION FUNCTION ##
#############################

def snpextraction(human1,mouse1,totalrows):
	for i,base in enumerate(human1): # REMOVE GAPS IN HUMAN - LENGTH EQUAL TO REGION WIDTH
		if base == mouse1[i]:
			pos = start + i
			continue
		if base == "-":
			continue
		if base != mouse1[i]:
			pos = start + i
			newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+mouse1[i]+'\t'+'.\t.\t.\t'+'GT'+'0'+'1'
			totalrows = totalrows + newrow	
	return(totalrows)

############################
## READING THE INPUT FILE ##
############################

with open(args.input,'r') as file:
	for line in file:
		if line[0]==">": 
			trigger+=1
		if trigger==3: # Every Human (or reference species) line after the first
		# REMOVE THE GAPS FROM THE HUMAN SEQUENCE:
			totalrows = snpextraction(human1,mouse1,totalrows)
			trigger=1 # Sets the trigger to 1 again					
			human1 = mouse1 = ""
			headers = []
			counter+=1
		# STORES THE POSITION OF THE CURRENT FRAGMENT:
		if trigger == 1 and line[0]==">": # HUMAN HEADER (including first)
			match = re.search(r'(chr\d+:)(\d+)-(\d+)',line)	# I've replaced \d\d with \d+
			if match: 
				start = int(match.group(2))
				end1 = int(match.group(3))
		# STORE INFORMATION TO CONSTRUCT THE FASTA FILE:	
			headers.append(line)
		elif trigger == 2 and line[0]==">": # MOUSE HEADER
			headers.append('\n'+line)
		elif trigger == 1 and line[0] != ">": # HUMAN SEQUENCE
			human1=human1+line.strip()
		elif trigger == 2 and line[0] != ">": # MOUSE SEQUENCE
			mouse1=mouse1+line.strip()

#################################
## CREATING THE FINAL VCF FILE ##
#################################

out = open(args.output,'w') # Open the final output file. 
out.write('##fileformat=VCFv4.1\n') # We assume it is always going to be the same version (but it can be modified)
out.write('##contig=<ID=%s,length=%d>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' % (args.chrom,seqlength))
out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHuman\tMouse')			