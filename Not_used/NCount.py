#!/usr/bin/python
seq = ""
import re
import sys
filename = sys.argv[1]
with open(filename) as file:
	for line in file:
		if line[0]==">":
			continue
		m = re.search(r'[ACGTacgt]',line) # Search tries all possible starting points in the string. Match just tries the first starting point.
		if m:
			print("Unique sequence: %s") %line # Adds any Ns before the match
			seq = seq + line[0:(m.start())]
			break
		else: 	
			seq = seq + line.rstrip()
print(len(seq)) # Amount of Ns in the sequence
print(seq[-10:])

			



