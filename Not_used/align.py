#!/usr/bin/env python3

# This scripts creates a set containing all the regions in the MFA alignment. Next, it subtracts this set from another set
# covering the entire sequence (start, end) to identify the unaligned regions.

import sys
import re
import gc

##########################
## VARIABLE DECLARATION ##
##########################

chrini = int()
chrfin = int()
aligned = set()
interval = set()

####################################
## DEFINTIION OF ALIGNED INTERVAL ##
####################################

with open(sys.argv[1]) as file:
	for line in file:	
		match = re.search(r'chr\d+:(\d+)-(\d+)',line)
		if match: 
			instart = int(match.group(1))
			inend = int(match.group(2))
		# IDENTIFY THE START AND THE END OF THE ENTIRE ALIGNED REGION:	
		if chrini == 0: 
			chrini = instart
		if inend > chrfin:
			chrfin = inend
		# ADD THE ALIGNED REGION TO THE TOTAL ALIGNMENT SET:
		aligned.update(set(range(instart,inend+1)))	# 1) Update for nclude the last position of the interval

########################################
## DEFINTIION OF NON-ALIGNED INTERVAL ##
########################################

unaligned = set(range(chrini,chrfin))-aligned
del aligned
gc.collect()
print("Region covered:",chrini,"-",chrfin)
print(sorted(unaligned))