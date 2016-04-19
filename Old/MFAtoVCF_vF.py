#!/usr/bin/env python3

# MFAtoVCF, Version 1.1B - This version replaces the 1 in the mouse column with dots when the alternative allele is missing (.).

# This script converts multi-FASTA alignments (.mfa) files structured as pairs of aligned sequences, each of which covers a region
# of a chromosome, into VCF files containing SNP data. For overlapping regions of the reference species that may match different parts
# of the target species genome, the following criteria are applied:
# 1) Any non-divergent position -> No SNP
# 2) Same divergence, a single SNP -> SNP
# 3) 2 different divergences-> Missing (.) 
# 4) Mouse gap -> Missing (.) 

# WARNING: This program must be provided with the reference sequence of that chromosome. 
# WARNING: By default, the script does NOT convert non-aligned regions to missing information. To do so, specif the "-n" option.
# UPDATES 1: The script has been cleaned of lines added during the testing stage.
# UPDATES 2: The "Human" column has been removed so that it does not interfere with the final merge.
# UPDATES 3: Now the VCF file is automatically converted to VCF.GZ with bgzip and indexed with tabix.
# UPDATES 4: The alternative allele is diploid for more uniformity with the 1000 GP (1/1 or ./.).
# FUTURE UPDATES: Study the feasibility of adding parallelism (multiprocessing module).

import time
import argparse
import re
import subprocess

t0 = time.clock()

#######################
## PARSING ARGUMENTS ##
#######################

parser = argparse.ArgumentParser(description='Extract the variants from a multi-fasta (MFA) alignment in VCF format. Any position absent from the alignment or ambiguous is considered missing (.).  By default, the tool compresses the VCF file with bgzip and indexes it with tabix, Currently it only handles single .mfa files containing information from one chromosome. ')
parser.add_argument('input',type=str,help="Input file, in MFA format")
parser.add_argument('-s','--sequence',type=str,required=True,help="DNA sequence of the chromosome")
parser.add_argument('-o','--output',type=str,default="outdef",help="Output file, in VCF format")
parser.add_argument('-c','--chrom',type=str,default="22",help="Chromosome identifier")
parser.add_argument('-np','--noprocess',action='store_true',help="Keeps the output in VCF format, skipping compression and indexing")
parser.add_argument('-t','--test',action='store_true',help="Testing mode, keeps temporal files for examination")

args = parser.parse_args()

if args.output == "outdef":
	args.output = args.input.rsplit('.')[0]+".vcf"

###########################
## IMPORT FASTA SEQUENCE ##
###########################

sequence = open(args.sequence, "rU") # r = read, U = Universal
refgenome = ""
upper = str.upper
for line in sequence: 
	if line[0] != ">":
		refgenome = refgenome + upper(line.rstrip())	

# Python has a module for importing FASTA sequences, but it is slower than this loop:
# from Bio import SeqIO
# refgenome = SeqIO.read(sequence,'fasta') # SeqIO.read expects only one FASTA record
# refgenome = refgenome.upper()

############################
## INITIALIZING VARIABLES ##
############################

# VCF MANIPULATION:
start0,end0 = 0,0
trigger = 0 # Trigger=0 - Start; Trigger=1 - 1st human seq; Trigger=2 - 1st mouse seq; Trigger=3 - 2nd human seq
arrayvcf = [None]*len(refgenome) # Must have the same size as the genome to compare each position
human,mouse = "",""

# VCF FILE HEADER:
headers = []
chrom = args.chrom
counter = 0 # Count the number of processed FASTA regions

#########################
## COMPARISON FUNCTION ##
#########################

def comparison(human,mouse,arrayvcf):
	pos = start-2 # Python is 0-based, but the reference genome is 1-based. Also, we add +1 at the beginning of the loop.
	for i,mnt in enumerate(mouse):
		if human[i] == "-":
			continue
		pos = pos + 1	
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
			arrayvcf = comparison(human,mouse,arrayvcf)
			trigger=1 # Sets the trigger to 1 again					
			human,mouse = "",""
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
			human=human+line.strip()
		elif trigger == 2 and line[0] != ">": # MOUSE SEQUENCE
			mouse=mouse+line.strip()
		if args.test and line[0] == ">":
			print (headers[0].strip())	
 	# UPON REACHING THE END OF THE FULE:
	arrayvcf = comparison(human,mouse,arrayvcf)
	human = mouse = ""
	headers = []
	counter+=1				

#################################
## CREATING THE FINAL VCF FILE ##
#################################

out = open(args.output,'w') # Open the final output file. 
out.write('##fileformat=VCFv4.1\n') # We assume it is always going to be the same version (but it can be modified)
out.write('##contig=<ID=%s,length=%d>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' % (args.chrom,end0-start0))
out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMouse')			
for ref,variant in enumerate(arrayvcf[start0:end0]):
	if variant == None: # Non-aligned region: missing mouse information.
		out.write('\n%s\t%d\t.\t%s\t.\t.\t.\t.\tGT\t./.' % (chrom,ref+start0+1,refgenome[ref+start0]))
	elif variant == "." and variant != None: # Missing mouse information: variant = .
		out.write('\n%s\t%d\t.\t%s\t%s\t.\t.\t.\tGT\t./.' % (chrom,ref+start0+1,refgenome[ref+start0],variant))
	elif variant != 0 and variant != None: # Divergence, not polymorphism (0)
		out.write('\n%s\t%d\t.\t%s\t%s\t.\t.\t.\tGT\t1/1' % (chrom,ref+start0+1,refgenome[ref+start0],variant))	
out.close()

###################
## FINAL TOUCHES ##
###################

if args.noprocess == False:

	subprocess.call(['bgzip', args.output]) # The -c argument keeps the original file and prints the new one to screen
	subprocess.call(['tabix', '-p', 'vcf',args.output+'.gz'])

print("\nExecution complete in %d second(s). %d FASTA pair(s) processed covering a total of %s nucleotides." % (time.clock()-t0,counter,end0-start0))

exit()
