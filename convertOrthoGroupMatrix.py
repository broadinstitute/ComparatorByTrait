import os
import sys
from collections import defaultdict

orthofinder_gene_count_matrix = sys.argv[1]
strain_counts = defaultdict(int)

samples = []
data = [] 
with open(orthofinder_gene_count_matrix) as of:
	for i, line in enumerate(of):
		line = line.strip()
		ls = line.split('\t')
		if i == 0:
			samples = [x.split('.annotation')[0] for x in ls[1:-1]]	
			data.append(['#OrthoFinder'] + samples)
		else:
			og = ls[0]
			for j, val in enumerate(ls[1:-1]):
				s = samples[j]
				strain_counts[s] += int(val)
			data.append(ls[:-1])
printlist = []	
for i, s in enumerate(samples):
	val = ''
	if i == 0:
		val = '#'
	val = val + s + '=' + str(strain_counts[s])
	printlist.append(val)

print('\t'.join(printlist))
for l in data:
	print('\t'.join(l))

