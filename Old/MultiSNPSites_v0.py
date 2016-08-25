#!usr/bin/python

# This script...
# ADDITIONAL DESCRIPTION: The tool SNP-sites has to be in the same folder

import argparse
import subprocess
import re

#######################
## PARSING ARGUMENTS ##
#######################

parser = argparse.ArgumentParser(description='Obtain a VCF file from a MFA using SNP-sites')
parser.add_argument('-i','--input',type=str,required=True,help="Input file, in MFA format")
parser.add_argument('-o','--output',type=str,default="output.vcf",help="Output file, in VCF format")
parser.add_argument('-s','--start',type=int,default=10000,help="Window size for metrics calculation (default 10000)")
parser.add_argument('-c','--chrom',default="22",help="Step size for metrics calculation (default 0)")

args = parser.parse_args()

#####################################
## CREATING VCF FILES BY SEQ PAIRS ##
#####################################

counter = 0
header = ""
temp = ""
total = 0
out = open(args.output,'w') # Open the final output file. 
with open(args.input,'r') as file: 
	for line in file:
		if line[0]==">": # Also: if "human" is in line
			counter+=1
			print (counter)
		if counter==3:
			tfile = open("tempin.fasta",'w') # Stores the pair of sequences in a temporal file
			tfile.write(temp) # Writes the temp variable in a .fasta file
			tfile.close()
			subprocess.call(['snp-sites','-v','-o','tempout.vcf','tempin.fasta']) # Calls snp-sites and generates a VCF output
			
			# MANIPULATE VCF FILE
			with open("tempout.vcf") as tempvcf: # Opens the VCF output
				for row in tempvcf:
					length = re.search(r'length=(\d+)>',row)
					if length: 
						print(length.group(1))
						total = total + int(length.group(1))
					if row[0]!= "#":
						rsplit = row.split()
						pos = int(header) + int(rsplit[1]) -1 # Arithmetical sum!!! We 
						newrow = args.chrom + '\t' + str(pos) + '\t' + '\t'.join(rsplit[2:]) + '\n'
						print(newrow)
						out.write(newrow) # Appends the temporal VCF to the final output
						counter=1 # Sets the counter to 1 again

			temp = "" # Resets the temp variable to empty
		if counter == 1:
			match = re.search(r'(chr\d\d?:)(\d+)-',line)	
			if match: header = match.group(2)
		temp=temp+line 		
