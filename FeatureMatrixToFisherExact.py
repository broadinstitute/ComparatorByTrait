import os
import sys
import argparse
from collections import defaultdict
import scipy.stats as stats
import numpy as np
import math
import uuid

try:
	assert (sys.version_info[0] == 3)
except:
	sys.stderr.write("Please use Python-3.4 to run this program. Exiting now ...\n");
	sys.exit(1)

### Hardcoded directory for R script to perform Adjusting using Storey's q-value method/function.

RSCRIPT_PATH=os.path.dirname(os.path.realpath(__file__)) + '/qvalue_computation.R'
def p_adjust_bh(p):
	"""
	Benjamini-Hochberg p-value correction for multiple hypothesis testing.
	"""
	p = np.asfarray(p)
	by_descend = p.argsort()[::-1]
	by_orig = by_descend.argsort()
	steps = float(len(p)) / np.arange(len(p), 0, -1)
	q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
	return q[by_orig]

def q_adjust(p):
	"""
	Compute q-values from list of p-values using FDR correction method by Storey et al 2001.
	"""
	tmp_pvalues = "tmp.pvalues." + uuid.uuid4().hex
	tmp_qvalues = "tmp.qvalues." + uuid.uuid4().hex
	tmp_pvalues_handle = open(tmp_pvalues, 'w')
	for i, pv in enumerate(p):
		tmp_pvalues_handle.write(str(i) + '\t'  +str(pv) + '\n')
	tmp_pvalues_handle.close()
	os.system('Rscript %s %s %s' % (RSCRIPT_PATH, tmp_pvalues, tmp_qvalues))
	qvalues = []
	for i, line in enumerate(open(tmp_qvalues)):
		if i == 0: continue
		line = line.rstrip('\n')
		ls = line.split()
		qvalues.append(float(ls[1]))
	os.system('rm %s %s' % (tmp_pvalues, tmp_qvalues))
	return qvalues

def computeLog2Ratio(Aratio, Bratio):
	log2ratio = None
	if Aratio == 0 and Bratio == 0:
		log2ratio = "NA"
	elif Aratio > 0 and Bratio == 0:
		log2ratio = "+INF"
	elif Aratio == 0 and Bratio > 0:
		log2ratio = "-INF"
	else:
		log2ratio = str(math.log2(Aratio/Bratio))
	return log2ratio

def runAnalysis(matrix, groupA, groupB, output, qvalue):
	try: assert(os.path.isfile(matrix) and os.path.isfile(groupA) and os.path.isfile(groupB))
	except:
		sys.stderr.write("Unable to verify one of the input files is actually an existant file within reach. Exiting now ...\n")
		raise RuntimeError

	groupAstrains = set([strain.strip() for strain in open(groupA).readlines()])
	groupBstrains = set([strain.strip() for strain in open(groupB).readlines()])

	#print(groupAstrains)
	#print(groupBstrains)
	try: assert(len(groupAstrains.intersection(groupBstrains)) == 0)
	except:
		sys.stderr.write("Strainlists for Group A and Group B overlap, this is unaccepted for a Fisher Exact Test. Exiting now ...\n")
		raise RuntimeError

	feature_data = {}
	features = []
	pvalues = []
	total_feats = {}
	col_to_samp = {}
	groupAtotalsum = 0
	groupBtotalsum = 0
	groupAheader = []
	groupBheader = []
	for i, line in enumerate(open(matrix)):
		line = line.rstrip('\n').strip()
		if i == 0:
			line = line[1:]
			ls = line.split('\t')
			for val in ls:
				vs = ['='.join(val.split('=')[:-1]), val.split('=')[-1]]
				total_feats[vs[0]] = int(vs[1])
				if vs[0] in groupAstrains: groupAtotalsum += int(vs[1])
				elif vs[0] in groupBstrains: groupBtotalsum += int(vs[1])
		elif i == 1:
			line = line[1:]
			ls = line.split('\t')
			for j, val in enumerate(ls[1:]):
				col_to_samp[j] = val
				if val in groupAstrains: groupAheader.append(val)
				elif val in groupBstrains: groupBheader.append(val)
		else:
			ls = line.split('\t')
			feature = ls[0]
			groupAcounts = []
			groupBcounts = []
			groupAsum = 0
			groupBsum = 0
			wAsum = 0.0
			wBsum = 0.0
			for j, val in enumerate(ls[1:]):
				samp = col_to_samp[j]
				val = float(val)
				if samp in groupAstrains: 
					groupAsum += int(val); 
					groupAcounts.append(val);
					if total_feats[samp] != 0: 
						wAsum += float(val)/total_feats[samp]
				elif samp in groupBstrains: 
					groupBsum += int(val); 
					groupBcounts.append(val); 
					if total_feats[samp] != 0:
						wBsum += float(val)/total_feats[samp]

			if groupAsum > 0 or groupBsum > 0:
				groupAnzcounts = [x for x in groupAcounts if x > 0]; groupBnzcounts = [x for x in groupBcounts if x > 0]
				A_nz_median = 0.0; B_nz_median =0.0; A_nz_Q1 = 0.0; B_nz_Q1 = 0.0
				if groupAsum >0 :
					A_nz_Q1 = np.percentile(groupAnzcounts, 25)
					A_nz_median = np.percentile(groupAnzcounts, 50)
				if groupBsum > 0:
					B_nz_median = np.percentile(groupBnzcounts, 50)
					B_nz_Q1 = np.percentile(groupBnzcounts, 25)
				A_proportion = float(len([x for x in groupAcounts if x >= A_nz_Q1 and x > 0]))/len(groupAcounts)
				B_proportion = float(len([x for x in groupBcounts if x >= B_nz_Q1 and x > 0]))/len(groupBcounts)
				if True:
					groupAcompsum = groupAtotalsum - groupAsum
					groupBcompsum = groupBtotalsum - groupBsum
					confusion_matrix = [[groupAsum, groupBsum], [groupAcompsum, groupBcompsum]]
					oddsratio, pvalue = stats.fisher_exact(confusion_matrix)
					rankingStat = abs(A_proportion - B_proportion)*100.0
					if A_proportion < B_proportion: rankingStat = -rankingStat
					printlist = [A_nz_median, B_nz_median, A_proportion, B_proportion]
					feature_data[feature] = printlist
					features.append(feature)
					pvalues.append(pvalue)

	adj_pvalues = []
	if len(pvalues) > 0:
		if qvalue: adj_pvalues = q_adjust(pvalues)
		else: adj_pvalues = p_adjust_bh(pvalues)

	outf = open(output, 'w')
	outf.write('\t'.join(['#feature', 'rawPvalue', 'Q-value', 'NZ-Median-A', 'NZ-Median-B', 'Proportion-A', 'Proportion-B'])  + "\n")

	for i, f in enumerate(features):
		pvalue = pvalues[i]
		adj_pvalue = adj_pvalues[i]
		printlist = [f, str(pvalue), str(adj_pvalue)]
		printlist += [str(x) for x in feature_data[f]]
		outf.write('\t'.join(printlist) + '\n')
	outf.close()

if __name__ == '__main__':
		# Parse arguments.

	parser = argparse.ArgumentParser(description="""
	This program will compute the adjusted p-value for performing a Fisher exact test for each annotation feature (e.g. KEGG Pathway)
	between two strainlists given an annotation feature matrix with the following header format:
	
	indicate the total number of genes (or features) corresponding to each of the column headers in the first line.

	#IN=5378        OUT=14735
	#domain IN      OUT
	
	The second line has the column names, which should match the first line.
	""")

	parser.add_argument('-i', '--matrix', type=str, help='Provide the input annotation feature matrix with the appropriate headings.', required=True)
	parser.add_argument('-a', '--groupA', type=str, help='Provide file to strainlist for Group A.', required=True)
	parser.add_argument('-b', '--groupB', type=str, help='Provide file to strainlist for Group B.', required=True)
	parser.add_argument('-o', '--output', type=str, help='Provide output file name.', required=True)
	parser.add_argument('-q', '--qvalue', action='store_true', help='Perform q-value FDR adjustment of p-values instead of Benjamini-Hochberg FDR adjustment. Please make sure R with the qvalue library is available in path.', required=False, default=False)

	args = parser.parse_args()
	runAnalysis(args.matrix, args.groupA, args.groupB, args.output, args.qvalue)
