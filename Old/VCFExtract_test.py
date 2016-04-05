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
parser.add_argument('-o','--output',type=str,default="outdef",help="Output file, in VCF format")
parser.add_argument('-c','--chrom',type=str,default="22",help="Chromosome identifier")
parser.add_argument('-t','--test',action='store_true',help="Testing mode, keeps temporal files for examination")

args = parser.parse_args()

if args.output == "outdef":
	args.output = args.input.rsplit('.')[0]+".vcf"

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

def snpextraction(human1,mouse1,totalrows,seqlength,human0,mouse0,start0,end0,last=0):	
	pos = 0
	### IF OVERLAP BETWEEN HUMAN REGIONS ###		
	if start1 < end0: 
		pos = start0
		# NON-OVERLAPPING REGION:
		for i,base in enumerate(human0): 
			if pos >= start1: # Pos can be larger if start1 < start0 (double overlap)
				break
			elif pos<start1 and base == mouse0[i]:
				pos = pos + 1 # The current base is pos, the next one is pos+1
			elif pos<start1 and mouse0[i] == '-':	
				newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+'.'+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1'
				totalrows = totalrows + newrow
				pos = pos + 1
			elif pos<start1 and base != mouse0[i]:
				newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+mouse0[i]+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1'
				totalrows = totalrows + newrow	
				pos = pos + 1 	
			
		# OVERLAPPING REGION:	
		ini_overlap = max(start0-start1,0) # In case a long sequence overlaps two+ other sequences. Gaps?!
		human1 = human1[ini_overlap:]
		mouse1 = mouse1[ini_overlap:]
		for x,base in enumerate(human1):
			if pos > end0: # BEGINNING OR END?
				break			
			if args.test:
				print("Start0",start0,"End0",end0,"Position",pos)
				print(base,mouse1[x],mouse0[pos-start0],x,pos-start0,"[human,mouse1,mouse0,x,pos-start0]")
				print('Human1',human1)
				print('Mouse1',mouse1)
			if base == "-":
				continue # Base NOT taken into account for the position
			if base == mouse1[x] or base == mouse0[pos-start0]: # 1) Any non-divergent position, we consider it non-divergent
				print("Non-divergent")
				pos = pos + 1
				continue
			elif mouse1[x] == "-" or mouse0[pos-start0] == "-": # 4) Gap in overlapping region -> Missing (.)
				print("Gap in mouse")
				newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+'.'+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1'
				totalrows = totalrows + newrow
			elif base != mouse1[x] and mouse1[x] == mouse0[pos-start0]: # 2) Same divergence, a single SNP
				print("Same divergence")
				newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+mouse0[pos-start0]+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1'
				totalrows = totalrows + newrow
			elif base != mouse1[x] and mouse1[x] != mouse0[pos-start0]: # 3) 2 different divergences -> Missing (.)
				print("2 different divergences")
				newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+'.'+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1'
				totalrows = totalrows + newrow
			pos = pos + 1 # The current base is pos, the next one is pos+1 (end of all)				
	
		# END OF THE OVERLAPPING REGION
		if end1 >= end0: # If the current seq is longer than the previous -> It becomes the new 0
			human1 = human1[x:] # It cannot go after pos > end0 because the end of the loop can be reached without this condition!
			mouse1 = mouse1[x:]
			print(human1)
			mouse0,human0 = "","" # We empty the current "mouse0" variable
			for y,base in enumerate(human1):	# We only store the non-overlapping part (right end)
				if base != "-":
					mouse0 = mouse0 + mouse1[y]
					human0 = human0 + base
			start0 = end0+1 # The new start is the end of the overlapping region (or end of the previous 0 sequence)
			end0 = end1	
			seqlength = seqlength + len(human0)		
		elif end0 >= end1: # If the previous sequence is longer than the current one
			human0 = human0[pos-start0:]
			mouse0 = mouse0[pos-start0:]	
			start0 = pos # The previous loop ends after adding +1 to pos, so the new start matches pos. E.g. if A goes from 1 to 100 and B goes from 10 to 50, the new start will be pos = 51				
			seqlength = seqlength + len(human0)

	### IF NOT OVERLAP BETWEEN HUMAN REGIONS ###		
	else:	
		# PRINT THE STORED SEQUENCE
		for i,base in enumerate(human0): 
			pos = start0 + i
			if base == mouse0[i]:
				continue
			if base != mouse0[i]:
				newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+mouse0[i]+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1'
				totalrows = totalrows + newrow
			if mouse0[i] == "-": # We do not count it as polymorphism, but as "missing information"
				newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+'.'+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1'
				totalrows = totalrows + newrow
		mouse0,human0 = "","" # We empty the current "mouse0" variable
		# REMOVE GAPS IN HUMAN AND STORE IT
		for i,base in enumerate(human1): 
			if base != "-":
				mouse0 = mouse0 + mouse1[i]
				human0 = human0 + base
		end0 = end1
		start0 = start1
		seqlength = seqlength + len(human0)

	### IF LAST REGION OF THE FILE ###	
	if last == 1:
		for i,base in enumerate(human0): # PRINT THE STORED SEQUENCE
			pos = start0 + i
			if base == mouse0[i] or mouse0[i] == "-":
				continue
			if base != mouse0[i]:
				newrow = '\n'+chrom+'\t'+str(pos)+'\t'+'.'+'\t'+base+'\t'+mouse0[i]+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1'
				totalrows = totalrows + newrow	
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
			human1,mouse1 = "",""
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
	totalrows,seqlength,human0,mouse0,start0,end0 = snpextraction(human1,mouse1,totalrows,seqlength,human0,mouse0,start0,end0,last=1)	
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

print("Execution complete. %d FASTA pair(s) processed covering a total of %s nucleotides." % (counter,seqlength))
exit()