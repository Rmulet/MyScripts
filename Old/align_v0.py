#!/usr/bin/env python3

import sys
import re

# Use sets instead of lists? For overlap between multiple regions of aligned/naligned
	
taligned = [] # Total aligned 
tnaligned = [] # Total no aligned
naligned = []
end0 = 0

with open(sys.argv[1]) as file:
	for line in file:	
		match = re.search(r'chr\d+:(\d+)-(\d+)',line)	# I've replaced \d\d with \d+
		if match: 
			start1 = int(match.group(1))
			end1 = int(match.group(2))
		aligned = list(range(start1,end1))
		if start1 not in taligned and taligned != []:
			naligned = list(range(end0,start1))
		taligned.append(aligned)
		tnaligned.append(naligned)
		aligned,naligned = [],[]
		start0,end0 = start1,end1
	print(tnaligned)	