import os
import sys
import argparse

try:
    assert (sys.version_info[0] == 3)
except:
    sys.stderr.write("Please use Python-3.4 to run this program. Exiting now ...\n");
    sys.exit(1)

def create_directory(directory):
    if not os.path.isdir(directory):
        os.system('mkdir %s' % directory)
    else:
        sys.stderr.write('Unable to create directory %s, it already exists. Overwriting!\n' % directory)

def create_parser():
    """ Parse arguments """
    parser = argparse.ArgumentParser(description="""
    Perform enrichment analysis for a given comparisons listing (see wiki page:).
    """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-c', '--comparisons', help='provide comparisons file.', required=True)
    parser.add_argument('-f', '--feature_matrix', help='feature matrix in appropriate format.', required=True)
    parser.add_argument('-o', '--output', help='output directory of results.', required=True)
    args = parser.parse_args()
    return args

# parse arguments
args = create_parser()

comparisons_file = args.comparisons
feature_matrix = args.feature_matrix
output_directory = os.path.abspath(args.output) + '/'

create_directory(output_directory)

strainlists_dir = output_directory + 'strainlists/'
create_directory(strainlists_dir)

enrichments_dir = output_directory + 'enrichments/'
create_directory(enrichments_dir)

scripts_directory = os.path.dirname(os.path.realpath(__file__))
enrichment_program = scripts_directory + '/FeatureMatrixToFisherExact.py'

tmp = {}
comp_id_to_name = {}
with open(comparisons_file) as ocf:
    for line in ocf:
        line = line.strip()
        if not line: continue
        ls = line.split('=')
        if line.startswith('//') and len(tmp) > 0:
            gAf = strainlists_dir + tmp['ID'] + '_groupA.txt'
            gBf = strainlists_dir + tmp['ID'] + '_groupB.txt'
            ogAf = open(gAf, 'w')
            ogBf = open(gBf, 'w')
            ogAf.write('\n'.join(tmp['GroupA'].split(',')) + '\n')
            ogBf.write('\n'.join(tmp['GroupB'].split(',')) + '\n')
            ogAf.close(); ogBf.close()
            enrichment_result = enrichments_dir + tmp['ID'] + '.txt'
            enrichment_cmd = ['python', enrichment_program, '-i', feature_matrix, '-a', gAf, '-b', gBf, '-o', enrichment_result]
            os.system(' '.join(enrichment_cmd))
            comp_id_to_name[tmp['ID']] = tmp['Name']
            tmp = {}
        else:
            tmp[ls[0]] = '='.join(ls[1:])

if len(tmp) > 0:
    gAf = strainlists_dir + tmp['ID'] + '_groupA.txt'
    gBf = strainlists_dir + tmp['ID'] + '_groupB.txt'
    ogAf = open(gAf, 'w')
    ogBf = open(gBf, 'w')
    ogAf.write('\n'.join(tmp['GroupA'].split(',')) + '\n')
    ogBf.write('\n'.join(tmp['GroupB'].split(',')) + '\n')
    ogAf.close();
    ogBf.close()
    enrichment_result = enrichments_dir + tmp['ID'] + '.txt'
    enrichment_cmd = ['python', enrichment_program, '-i', feature_matrix, '-a', gAf, '-b', gBf, '-o', enrichment_result]
    os.system(' '.join(enrichment_cmd))
    comp_id_to_name[tmp['ID']] = tmp['Name']

final_results = open(output_directory + 'consolidated_results.txt', 'w')
final_filtered_results = open(output_directory + 'consolidated_results.filt.txt', 'w')
for j, f in enumerate(os.listdir(enrichments_dir)):
    comparison_name = comp_id_to_name[f.split('.txt')[0]]
    with open(enrichments_dir + f) as of:
        for i, line in enumerate(of):
            if i == 0 and j == 0:
                final_results.write('comparison\t' + line)
                final_filtered_results.write('comparison\t' + line)
                continue
            elif i == 0: continue
            line = line.strip()
            ls = line.split('\t')
            qvalue = float(ls[2])
            groupA_prop = float(ls[-2])
            groupB_prop = float(ls[-1])
            final_results.write(comparison_name + '\t' + line + '\n')
            if qvalue < 0.05 and ( (groupA_prop >= 0.75 and groupB_prop <= 0.25) or (groupA_prop <= 0.25 and groupB_prop >= 0.75) ):
                final_filtered_results.write(comparison_name + '\t' + line + '\n')
final_results.close()
final_filtered_results.close()
