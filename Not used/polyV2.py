#!/usr/bin/python

# IMPORTANT: THIS CODE HAS BEEN WRITTEN IN PYTHON 3. IT WILL NOT WORK PROPERLY ON PYTHON 2. 

# This script calculates 3 nucleotide variation metrics: number of segregating sites, Watterson theta estimator and nucleotide diversity (Pi):
# (1) Reads the sequences in the 'input' file and stores them in the 'seqarray' list. Then it measures m and n on this array.
# (2) Counts the number of nucleotides of each class at a given position across all the sequences. That is, the number of A, C, G and T.
# (3) Calculates heterozigosity based on the count.
# i) If all nt are of one class, the position is homozygous for the sample analyized
# ii) If the sum of ACGT < n, some nucleotides are not A/C/G/T and the position is skipped
# iii) k, or number of differences, can be calculated as a sum of products of A,C,G,T --> stored in a new list
# (4) Calculates nucleotide variation metrics, defined as independent functions.
# (5) Prints the result to an output file in GFF format (with the aid of tabulators for each field)

import math
import argparse

#######################
## PARSING ARGUMENTS ##
#######################

parser = argparse.ArgumentParser(description='Calculate nucleotide variation metrics')

parser.add_argument('-i','--input',required=True,help="Input file, in FASTA format")
parser.add_argument('-o','--output',default="output.gff",help="Output file, in GFF format")
parser.add_argument('-w','--window',type=int,default=10000,help="Window size for metrics calculation (default 10000)")
parser.add_argument('-s','--step',default="window",help="Step size for metrics calculation (default 0)")
parser.add_argument('-r','--range',nargs=2,default=[10000001,11000000],type=int,help="Starting and end position in the genome (default 10000001-11000000)")
parser.add_argument('-d','--id',default="2L",type=str,help="ID of the sequence used to establish the coordinate system (e.g. chromosome)")

args = parser.parse_args()

fname = args.input
initial = args.range[0]
final = args.range[1]
window = args.window
if args.step == "window": # By default, the step is of the same size as the window
	step = window
else:
	step = int(args.step) 
seqid = args.id #  Used to generate the 'output' file
	
########################################
## READING FILE AND MEASURING N AND M ##
########################################
	
with open(fname) as file:
	seqarray = []
	for line in file:
		if line [0] != '>':
			seqarray.append(line.rstrip())
			
mtotal = len(seqarray[1]) # Number of nucleotides analyzed
n = len(seqarray) # Sample size (number of sequences). Invariable.

##################################################
## COUNTING THE NUMBER OF DIFFERENT NUCLEOTIDES ##
##################################################

# We create a vector: nt = (a,c,g,t) for each position. Then we add one to the corresponding nucleotide every time it is found; e.g. if A, nt[0]++1.
# This computes the number of 
# variations in a given position. If the result is 0, then k is also 0 (no differences)

totcount = []
for pos in range (mtotal):
	count = [0,0,0,0]
	for seq in seqarray:
		if (seq[pos]=='A' or seq[pos]=='a'):
			count[0] += 1
		elif (seq[pos]=='C' or seq[pos]=='c'):
			count[1] += 1
		elif seq[pos]=='G' or seq[pos]=='g':
			count[2] += 1
		elif seq[pos]=='T' or seq[pos]=='t':
			count[3] += 1
		else: # Therefore, in positions where there are no nt, the sum of count is lower than n
			break
	totcount.append(count)

##################################	
## NUCLEOTIDE VARIATION METRICS	##
##################################		
	
# Number of segregating sites per nucleotide:
	
def Number(m,k):
	S = m - k.count(0) # Number of segregating sites: total sites minus those that have no variation
	Sm = S/m
	return(Sm)

# Watterson theta estimator:

def Theta(m,k):
	S = m - k.count(0) # Number of segregating sites: total sites minus those that have no variation
	Sm = S/m
	sumat = 0
	for i in range (1,n): # Does not include the last number of the range (n)
		sumat = sumat + 1/i
	theta = Sm/sumat
	return(theta)

# Expected nucleotide heterozygosity (Pi)

def Pi(m,k):
	ktotal = sum(k)
	comb = math.factorial(n)/(2*math.factorial(n-2))
	pi = ktotal/(m*comb)
	return(pi)
			
#############################
## SLIDING WINDOW ANALYSIS ##
#############################	

numb_list = []; theta_list = []; pi_list = []
stwindows = [initial] # Vector that contains the starting position of every window

# Ktotals contain the total of variations detected in 

for i in range (1,mtotal,step):
	subcount = totcount[i-1:i+window-1] # We subtract 1 because Python is 0-based
	stwindows.append(stwindows[-1]+step) # Each window is separated by a distance = step
	msub = window
	ksub = [] 
	for nt in subcount: # Msub is m (number of nucleotides analyzed) in that window
		if sum(nt) < n:
			msub = msub-1
			continue
		nt.sort(reverse=True)
		k = nt[0]*sum(nt[1:3])+nt[1]*sum(nt[2:3])+nt[2]*sum(nt[3:3]) # Computes the number of heterozygous out of all pairs:
		ksub.append(k)
	
	if msub == 0: # In regions with large stretches of "N" or non-ACGT nucleotides, msub becomes 0 and metrics cannot be calculated (NA)
		numb_list.append(".")
		theta_list.append(".")
		pi_list.append(".")
	else: # If metrics can be calculated, they are stored in lists to be printed in the GFF file
		numb_list.append(Number(msub,ksub))
		theta_list.append(Theta(msub,ksub))
		pi_list.append(Pi(msub,ksub))

##################################
## PRINTING THE RESULTS TO FILE ##
##################################

out = open(args.output,'w') # Open the output file. By default, args.output = "output.gff"
	
out.write("##gff-version 3\n") 
for i,value in enumerate(numb_list):
	out.write("%s\tS\tbin\t%s\t%s\t%s\t.\t.\tName=Segregating sites per nucleotide\n" % (seqid,stwindows[i],stwindows[i]+window-1,value))
for i,value in enumerate(theta_list):
	out.write("%s\tTheta\tbin\t%s\t%s\t%s\t.\t.\tName=Waterson theta estimator\n" % (seqid,stwindows[i],stwindows[i]+window-1,value))	
for i,value in enumerate(pi_list):
	out.write("%s\tPi\tbin\t%s\t%s\t%s\t.\t.\tName=Nucleotide diversity\n" % (seqid,stwindows[i],stwindows[i]+window-1,value))	
	
out.close()