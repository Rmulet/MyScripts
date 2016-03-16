#!usr/bin/python

import re
import argparse 

#######################
## PARSING ARGUMENTS ##
#######################

parser = argparse.ArgumentParser(description='Obtain a VCF file from a MFA using the SNP-sites tool. Currently it is only prepared to handle single .mfa files containing information from one chromosome.')
parser.add_argument('-i','--input',type=str,required=True,help="Input file, in MFA format")
parser.add_argument('-p','--position',type=int,nargs='+',help="Starting position(s) of the region(s) to be extracted (no separators)") 
parser.add_argument('-t','--target',type=str,default="Human",help="Common name of the target species (check the file)")
parser.add_argument('-q','--query',type=str,default="Mouse",help="Common name of the target species (check the file)")
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
		match1 = re.search(r'%s.*chr\d+:%s'% (args.target,args.position[order]),line) # %s 
		match0 = re.search(r'%s.*chr\d+:\d+-' % args.query,line)
		if match1:
			slicer = 1
			if len(args.position) > 1 and len(args.position) > order+1: # If we define more than one position
				order+=1
		elif match0:
			slicer = 0
		if slicer == 1:
			print (line)	