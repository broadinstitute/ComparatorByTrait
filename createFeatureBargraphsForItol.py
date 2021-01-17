import os
import sys
from collections import defaultdict
import math
import argparse

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Create 
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-f', '--feature_matrix', help='Feature matrix', required=True)
	parser.add_argument('-o', '--output', help='output_directory', required=True)
	args = parser.parse_args()
	return args

myargs = create_parser()
feature_matrix = myargs.feature_matrix
resdir = os.path.abspath(myargs.output) + '/'

if not os.path.isdir(resdir):
	os.system('mkdir %s' % resdir)
else:
	sys.stderr.write('Output directory already exists! Exiting now ...\n')
	sys.exit(1)


feature_strain_dict = defaultdict(dict)
col_to_strain = {}
all_strains = set([])
all_features = set([])
for i, line in enumerate(open(feature_matrix)):
	line = line.rstrip('\n')
	ls = line.split('\t')
	if i == 0: continue
	elif i == 1:
		for j, val in enumerate(ls[1:]):
			col_to_strain[j] = val
			all_strains.add(val)
	else:
		kp = ls[0]
		all_features.add(kp)
		for j, val in enumerate(ls[1:]):
			feature_strain_dict[kp][col_to_strain[j]] = int(val)
		
for p in all_features:
	outf = open(resdir + p.replace('/', '') + '.txt', 'w')	
	outf.write('DATASET_MULTIBAR\nSEPARATOR TAB\nDATASET_LABEL\t'+ p + '\nCOLOR\t#ff0000\n')
	outf.write('FIELD_COLORS\t#4682B4\n')
	outf.write('FIELD_LABELS\t' +  p + '\n')
	outf.write('ALIGN_FIELDS\t1\nBORDER_WIDTH\t20\nBORDER_COLOR\t#ffffff\nDATA\n')
	for s in all_strains:
		pathway_val = str(feature_strain_dict[p][s])
		outf.write(s + '\t' + pathway_val + '\n')
	outf.close()