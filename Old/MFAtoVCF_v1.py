#!/usr/bin/env python3

# MFAtoVCF, Version 0.6

# This script converts multi-FASTA alignments (.mfa) files structured as pairs of aligned sequences, each of which covers a region
# of a chromosome, into VCF files containing SNP data. For overlapping regions of the reference species that may match different parts
# of the target species genome, the following criteria are applied:
# 1) Any non-divergent position -> No SNP
# 2) Same divergence, a single SNP -> SNP
# 3) 2 different divergences-> Missing (.) 
# 4) Mouse gap -> Missing (.) 

# WARNING: This program must be provided with the reference sequence of that chromosome. 
# UPDATES: The script has been optimized by removing the 'remgaps' function and integrating it in the main 'comparison' function. This 
# reduces the memory usage (there is only one pair of human-mouse variables) and increases speed.

import argparse
import re
import os

#######################
## PARSING ARGUMENTS ##
#######################

parser = argparse.ArgumentParser(description='Obtain a VCF file from a MFA using the SNP-sites tool. Currently it is only prepared to handle single .mfa files containing information from one chromosome.')
parser.add_argument('input',type=str,help="Input file, in MFA format")
parser.add_argument('-s','--sequence',type=str,required=True,help="DNA sequence of the chromosome")
parser.add_argument('-o','--output',type=str,default="outdef",help="Output file, in VCF format")
parser.add_argument('-c','--chrom',type=str,default="22",help="Chromosome identifier")
parser.add_argument('-n','--nonaligned',action='store_true',help="Non-aligned regions are printed as missing information")
parser.add_argument('-t','--test',action='store_true',help="Testing mode, keeps temporal files for examination")

args = parser.parse_args()

if args.output == "outdef":
	args.output = args.input.rsplit('.')[0]+"son.vcf"

###########################
## IMPORT FASTA SEQUENCE ##
###########################

from Bio import SeqIO
sequence = open(args.sequence, "rU") # r = read, U = Universal
refgenome = SeqIO.read(sequence,'fasta') # SeqIO.read expects only one FASTA record
refgenome = refgenome.upper()
# print(refgenome[50957596-1])

############################
## INITIALIZING VARIABLES ##
############################

# VCF MANIPULATION:
start0,end0 = 0,0
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

def comparison(human1,mouse1,arrayvcf):
	pos = start-2 # Python is 0-based, but the reference genome is 1-based. Also, we add +1 at the beginning of the loop.
	#print("start",start,"end",end)
	#print("Region",human2[1])
	#print("Refgenome",refgenome[pos:end-1])
	for i,mnt in enumerate(mouse1):
		if human1[i] == "-":
			continue
		pos = pos + 1	
		#print("Region",mnt)
		#print("Refgenome",refgenome[pos],pos)
		# FILLED POSITION:
		if arrayvcf[pos] != None:
			if mnt == arrayvcf[pos]: # 2) Same divergence, a single SNP
				continue
			elif mnt != arrayvcf[pos] and arrayvcf[pos] == 0: # 1) Any non-divergent position -> 0 (former non-divergent)
				continue
			elif mnt != arrayvcf[pos] and mnt == refgenome[pos]: # 1) Any non-divergent position -> 0 (new non-divergent)
				arrayvcf[pos] = 0	
			elif mnt != arrayvcf[pos] and mnt != refgenome[pos] and arrayvcf[pos] != 0: # 3) 2 different divergences & 4) Mouse gap -> Missing (.) 
				arrayvcf[pos] = "."
		# EMPTY POSITION:
		elif arrayvcf[pos] == None:
			if mnt == refgenome[pos]: # No divergence
				arrayvcf[pos] = 0
			elif mnt != refgenome[pos] and mnt != "-": # Divergence
				arrayvcf[pos] = mnt
			elif mnt != refgenome[pos] and mnt == "-": # Mouse gap
				arrayvcf[pos] = "."
	#print("arrayvcf",start,end,arrayvcf[start-1:end])
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
		if end > end0:
			end0 = end
		# STORE INFORMATION TO CONSTRUCT THE FASTA FILE:	
		elif trigger == 2 and line[0]==">": # MOUSE HEADER
			headers.append('\n'+line)
		elif trigger == 1 and line[0] != ">": # HUMAN SEQUENCE
			human1=human1+line.strip()
		elif trigger == 2 and line[0] != ">": # MOUSE SEQUENCE
			mouse1=mouse1+line.strip()
		if args.test and line[0] == ">":
			print (headers[0].strip())	
 	# UPON REACHING THE END OF THE FULE:
	arrayvcf = comparison(human1,mouse1,arrayvcf)
	human1 = mouse1 = ""
	headers = []
	counter+=1				

#################################
## CREATING THE FINAL VCF FILE ##
#################################

out = open(args.output,'w') # Open the final output file. 
out.write('##fileformat=VCFv4.1\n') # We assume it is always going to be the same version (but it can be modified)
out.write('##contig=<ID=%s,length=%d>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' % (args.chrom,100))
out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHuman\tMouse')			
for ref,variant in enumerate(arrayvcf[start0:end0]):
	if variant == None and args.nonaligned:
		out.write('\n'+chrom+'\t'+str(ref+start0+1)+'\t'+'.'+'\t'+refgenome[ref+start0]+'\t'+'N'+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1')
	if variant != 0 and variant != None:
		out.write('\n'+chrom+'\t'+str(ref+start0+1)+'\t'+'.'+'\t'+refgenome[ref+start0]+'\t'+variant+'\t'+'.\t.\t.\t'+'GT'+'\t'+'0'+'\t'+'1')
out.close()

print("Execution complete. %d FASTA pair(s) processed covering a total of %s nucleotides." % (counter,end0-start0))
exit()