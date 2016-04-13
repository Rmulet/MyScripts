#!/usr/bin/env python3

# This script extracts the nucleotide in a given position of a FASTA file or the sequence in a range of positions.
# Also, if the "-l" or "line" option is specified, it returns the nucleotide or nucleotides in the line or range of lines specified.
# WARNING: Ranges of positions must be indicated with either "-" or ",". Spaces are not accepted.
# UPDATE: The "-t" or "-tellme" flag reveals the line where a given nucleotide is found in the FASTA file.

import sys
import re
import math

def display_usage():
	print("\nThis script extracts the nucleotide in a given position of a FASTA file or the sequence in a range of positions.")
	print ("\nUsage:\n%s [FASTA sequence] [Options] [Positions] \n" % sys.argv[0])

# If less than two arguments supplied, display usage 
if len(sys.argv) < 2:
	print ("\nERROR: Missing arguments")
	display_usage()
	exit()

# Check whether user had supplied -h or --help . If yes display usage 
if sys.argv[1] == "-h" or sys.argv[1] == "--help":
	display_usage()
	exit()

sequence = open(sys.argv[1], "rU") # r = read, U = Universal
refgenome = ""
upper = str.upper
for line in sequence: 
	if line[0] != ">":
		refgenome = refgenome + upper(line.rstrip())	

# TO DETERMINE THE NUCLEOTIDE BASED ON THE LINE NUMBER:
if sys.argv[2] == "-l" or sys.argv[2] == "line":
	coords = re.split("-|,",sys.argv[3])
	if len(coords) == 2:
		pos = [(int(coords[0])-1)*50,(int(coords[1])-1)*50]# We subtract the header line
		print(pos)
		print(refgenome[pos[0]-1:pos[1]-1])  # -1 because Python is 0-based.
	else:
		pos = (int(sys.argv[3])-1)*50 # We subtract the header line
		print(pos)
		print(refgenome[pos-1])	# -1 because Python is 0-based

if sys.argv[2] == "-t" or sys.argv[2] == "tellme":
	line = math.ceil(int(sys.argv[3])/50)+1
	print(line)

# RETURN THE NUCLEOTIDE IN A POSITION OR RANGE THEREOF
else:
	coords = re.split(",|-",sys.argv[2])
	if len(coords) == 2:
		pos = [int(coords[0]),int(coords[1])]
		print(pos)
		print(refgenome[pos[0]-1:pos[1]-1]) # -1 because Python is 0-based
	else:
		pos = int(sys.argv[2])
		print(refgenome[pos-1]) # -1 because Python is 0-based

