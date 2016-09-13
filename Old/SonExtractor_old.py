#!/usr/bin/env python3

# SonExtracter, Version 0.0

import argparse
import re
import os

#######################
## PARSING ARGUMENTS ##
#######################

parser = argparse.ArgumentParser(description='Obtain a VCF file from a MFA using the SNP-sites tool. Currently it is only prepared to handle single .mfa files containing information from one chromosome.')
parser.add_argument('-i','--input',type=str,required=True,help="Input file, in MFA format")
parser.add_argument('-s','--sequence',type=str,required=True,help="DNA sequence of the chromosome")
parser.add_argument('-o','--output',type=str,default="outdef",help="Output file, in VCF format")
parser.add_argument('-c','--chrom',type=str,default="22",help="Chromosome identifier")
parser.add_argument('-t','--test',action='store_true',help="Testing mode, keeps temporal files for examination")

args = parser.parse_args()

if args.output == "outdef":
	args.output = args.input.rsplit('.')[0]+"2.vcf"

###########################
## IMPORT FASTA SEQUENCE ##
###########################

from Bio import SeqIO
sequence = open(args.sequence, "rU") # r = read, U = Universal
refgenome = SeqIO.read(sequence,'fasta') # SeqIO.read expects only one FASTA record
# print(refgenome[50957596-1])

############################
## INITIALIZING VARIABLES ##
############################

# VCF MANIPULATION:
start0 = 0
trigger = 0 # Trigger=0 - Start; Trigger=1 - 1st human seq; Trigger=2 - 1st mouse seq; Trigger=3 - 2nd human seq
arrayvcf = [None]*len(refgenome) # Must have the same size as the genome to compare each position
human1,mouse1 = "",""

# VCF FILE HEADER:
headers = []
chrom = args.chrom
counter = 0 # Count the number of processed FASTA regions

#########################
## COMPARISON FUNCTION ##
#########################

def remgaps(human1,mouse1):
	human2,mouse2 = "",""
	for i,base in enumerate(human1): # REMOVE GAPS IN HUMAN SEQUENCE
		if base != "-":
			mouse2 = mouse2 + mouse1[i]
			human2 = human2 + base
	return(human2,mouse2)		

def comparison(human1,mouse1,arrayvcf):
	pos = start-1 # Python is 0-based, but the reference genome is 1-based
	human2, mouse2 = remgaps(human1,mouse1)
	print("start",start,"end",end)
	#print("Region",human2[1])
	#print("Refgenome",refgenome[pos:end-1])
	for mnt in mouse2:
		print("Region",mnt)
		print("Refgenome",refgenome[pos])
		# FILLED POSITION:
		if arrayvcf[pos] != None:
			if mnt == arrayvcf[pos]: # 2) Same divergence, a single SNP
				continue
			elif mnt != arrayvcf[pos] and mnt != refgenome[pos]: # 3) 2 different divergences & 4) Mouse gap -> Missing (.) 
				arrayvcf[pos] = "."
			elif mnt != arrayvcf[pos] and mnt == refgenome[pos]: # 1) Any non-divergent position, we consider it non-divergent
				arrayvcf[pos] = refgenome[pos]
		# EMPTY POSITION:
		elif arrayvcf[pos] == None:
			if mnt == refgenome[pos]:
				arrayvcf[pos] = refgenome[pos]
			elif mnt != refgenome[pos]: 
				arrayvcf[pos] = mnt	
		pos = pos + 1
	print("arrayvcf",arrayvcf[start-1:end])
	return (arrayvcf)

############################
## READING THE INPUT FILE ##
############################

with open(args.input,'r') as file:
	for line in file:
		if line[0]==">": 
			trigger+=1
		if trigger==3: # Every Human (or reference species) line after the first
		# REMOVE THE GAPS FROM THE HUMAN SEQUENCE:
			arrayvcf = comparison(human1,mouse1,arrayvcf)
			trigger=1 # Sets the trigger to 1 again					
			human1,mouse1 = "",""
			headers = []
			counter+=1
		# STORES THE POSITION OF THE CURRENT FRAGMENT:
		if trigger == 1 and line[0]==">": # HUMAN HEADER (including first)
			match = re.search(r'(chr\d+:)(\d+)-(\d+)',line)	# I've replaced \d\d with \d+
			if match: 
				start = int(match.group(2))
				end = int(match.group(3))
		headers.append(line)
		if start0 == 0:
			start0 = start				
		# STORE INFORMATION TO CONSTRUCT THE FASTA FILE:	
		elif trigger == 2 and line[0]==">": # MOUSE HEADER
			headers.append('\n'+line)
		elif trigger == 1 and line[0] != ">": # HUMAN SEQUENCE
			human1=human1+line.strip()
		elif trigger == 2 and line[0] != ">": # MOUSE SEQUENCE
			mouse1=mouse1+line.strip()
 	# UPON REACHING THE END OF THE FULE:
	arrayvcf = comparison(human1,mouse1,arrayvcf)
	human1 = mouse1 = ""
	headers = []
	counter+=1				