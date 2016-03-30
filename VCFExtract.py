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
human2,mouse2 = "",""
totalrows = ""
human0 = ""
mouse0 = ""
start0,end0 = 0,0 # Longest end so far

# VCF FILE HEADER:
headers = [] # First human, second mouse
seqlength = 0
chrom = args.chrom
counter = 0 # Count the number of processed FASTA regions

#############################
## SNP EXTRACTION FUNCTION ##
#############################

def snpextraction(human1,mouse1,totalrows,seqlength,human0,mouse0,start0,end0):	
	pos = 0
	### IF OVERLAP BETWEEN HUMAN REGIONS ###		
	if start1 < end0: 
		pos = start0
		# NON-OVERLAPPING REGION:
		for i,base in enumerate(human0): 
			if (pos<start1 and base == human0[i]) or (pos<start1 and mouse0[i] == '-'):
				continue
			if 	pos<start1 and base != mouse0[i]:
				newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+mouse0[i]+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1'
				totalrows = totalrows + newrow		
			if pos == start1:
				break
			pos = pos + 1 # The current base is pos, the next one is pos+1
		# OVERLAPPING REGION:	
		for x,base in enumerate(human1):
			if base == "-":
				continue # Not taken into account
			if base == mouse1[x] or base == mouse0[pos-start0]: # 1) Any non-divergent position, we consider it non-divergent
				pos = pos + 1
				continue
			elif base != mouse1[x] and mouse1[x] == mouse0[pos-start0]: # 2) Same divergence, a single SNP
				newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+mouse0[i]+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1'
				totalrows = totalrows + newrow
			elif mouse1[x] == "-" or mouse0[pos-start0] == "-": # 4) Gap in overlapping region -> Missing (.)
				newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+'.'+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1'
				totalrows = totalrows + newrow
			elif base != mouse1[x] and mouse1[x] != mouse0[pos-start0]: # 3) 2 different divergences -> Missing (.)
				newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+'.'+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1'
				totalrows = totalrows + newrow
			if pos == end0:
				break
			pos = pos + 1 # The current base is pos, the next one is pos+1				
		# END OF THE OVERLAPPING REGION
		if end1 >= end0:	
			mouse0 = "" # We empty the current "mouse0" variable
			for base in enumerate(human1[pos-start1:]):	# We only store the non-overlapping part (right end)
				mouse0 = mouse0 + mouse1[pos-start1]
				human0 = human0 + base
				pos = pos + 1
			end0 = end1
			start0 = pos-start1
			seqlength = seqlength + len(human0)
	### IF NOT OVERLAP BETWEEN HUMAN REGIONS ###		
	else:			
		for i,base in enumerate(human0): # PRINT THE STORED SEQUENCE
			pos = start0 + i
			if base == mouse0[i] or mouse0[i] == "-":
				continue
			if base != mouse0[i]:
				newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+mouse0[i]+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1'
				totalrows = totalrows + newrow
		for i,base in enumerate(human1): # REMOVE GAPS IN HUMAN AND STORE IT
			mouse = "" # We empty the current "mouse0" variable
				if base != "-":
					mouse0 = mouse0 + mouse1[i]
					human0 = human0 + base
				end0 = end1
				start0 = pos-start1 
			seqlength = seqlength + len(human0)		
	return(totalrows,seqlength,human0,mouse0,start0,end0)

############################
## READING THE INPUT FILE ##
############################

with open(args.input,'r') as file:
	for line in file:
		if line[0]==">": 
			trigger+=1
		if trigger==3: # Every Human (or reference species) line after the first
		# REMOVE THE GAPS FROM THE HUMAN SEQUENCE:
			totalrows,seqlength,human0,mouse0,start0,end0 = snpextraction(human1,mouse1,totalrows,seqlength,human0,mouse0,start0,end0)
			trigger=1 # Sets the trigger to 1 again					
			human1 = mouse1 = ""
			if args.test:
				print (headers[0].strip()+':',len(human2))
			headers = []
			counter+=1
		# STORES THE POSITION OF THE CURRENT FRAGMENT:
		if trigger == 1 and line[0]==">": # HUMAN HEADER (including first)
			match = re.search(r'(chr\d+:)(\d+)-(\d+)',line)	# I've replaced \d\d with \d+
			if match: 
				start1 = int(match.group(2))
				end1 = int(match.group(3))
		# STORE INFORMATION TO CONSTRUCT THE FASTA FILE:	
			headers.append(line)
		elif trigger == 2 and line[0]==">": # MOUSE HEADER
			headers.append('\n'+line)
		elif trigger == 1 and line[0] != ">": # HUMAN SEQUENCE
			human1=human1+line.strip()
		elif trigger == 2 and line[0] != ">": # MOUSE SEQUENCE
			mouse1=mouse1+line.strip()
 	# UPON REACHING THE END OF THE FULE:
	totalrows,seqlength,human0,mouse0,start0,end0 = snpextraction(human1,mouse1,totalrows,seqlength,human0,mouse0,start0,end0)	
	if args.test:
		print (headers[0].strip()+':',len(human2))
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