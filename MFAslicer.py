#!/usr/bin/python3

import re
import argparse 

#######################
## PARSING ARGUMENTS ##
#######################

parser = argparse.ArgumentParser(description='Obtain a VCF file from a MFA using the SNP-sites tool. Currently it is only prepared to handle single .mfa files containing information from one chromosome.')
parser.add_argument('-i','--input',type=str,required=True,help="Input file, in MFA format")
parser.add_argument('-p','--position',type=int,nargs='+',help="Starting position(s) of the region(s) to be extracted (no separators)") 
parser.add_argument('-s','--species',type=str,default="Human",help="Common name of the target species (check the file)")
# nargs = '+' indicates that at least one argument is required, and they will be stored in a list.

args = parser.parse_args()

# REMINDER: Separate lists with commas and hyphens

############################
## INITIALIZING VARIABLES ##
############################

slicer = 0
order = 0

#########################
## EXTRACTING SEQUENCE ##
#########################

with open(args.input,'r') as file: 
	for line in file:
		match1 = re.search(r'%s.*chr\d+:%s-'% (args.species,args.position[order]),line) # 
		match0 = re.search(r'%s.*chr\d+:\d+-' % args.species,line)
		if match1:
			slicer = 1
			if len(args.position) > 1 and len(args.position) > order+1: # If we define more than one position
				order+=1
		if match0 and not match1:
			slicer = 0
		if slicer == 1:
			print (line.rstrip())	