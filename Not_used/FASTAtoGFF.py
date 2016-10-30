#!/usr/bin/python3

import re
import argparse 
import time
import os

#######################
## PARSING ARGUMENTS ##
#######################

parser = argparse.ArgumentParser(description='Replace positions in a FASTA sequence with a code corresponding to a GFF feature')
parser.add_argument('input',type=str,help="Input file, in txt format, containing the gene list")
parser.add_argument('-s','--sequence',type=str,required=True,help="DNA sequence of the chromosome")

args = parser.parse_args()

# REMINDER: Separate lists with commas and hyphens

###########################
## IMPORT FASTA SEQUENCE ##
###########################

'''sequence = open(args.sequence, "rU") # r = read, U = Universal # 0.81 seconds
refgenome = ""
upper = str.upper
for line in sequence: 
	if line[0] != ">":
		refgenome = refgenome + upper(line.rstrip())	'''
# gffseq = [0]*len(refgenome)


# We want the sequence! If not, we can just execute the following:
# os.system("grep -v '>' chr22.fa | tr '\n' ' ' | sed 's/\W//g' | wc -c") # 0.00031
ref=os.popen("grep -v '>' chr22.fa | tr '\n' ' ' | sed 's/\W//g' | wc -c").read().split() # Returns characters
gffseq = [0]*int(ref[0])


###########################
## MANIPULATING SEQUENCE ##
###########################

with open(args.input,'r') as file: 
	for line in file:
		coords = [int(n) for n in line.split()] # Extracts coordinates and converts into integer
		stretch = coords[1]-coords[0]+1
		gffseq[coords[0]:coords[1]+1] = stretch*[0]


###################
## FINAL TOUCHES ##
###################

out = open("output.txt",'w') # Open the final output file. 
out.write(''.join(map(str,gffseq)))
out.close()

#awk '$3 == "CDS"' chr22.gff		