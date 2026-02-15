# Enrichment analysis for HCP results

import pickle
import csv
import sys, os

from collections import defaultdict
import numpy as np
from localgraph import node_cluster, prune_graph

from enrichment_helpers import is_cortical_dk_feature, build_roi_to_system, yeo_enrichment

sys.path.insert(0, os.path.abspath('../../..'))
from applications.hcp.data.format_feature_names import format_feature_names
from methods import method_names, silggm_methods

#----------------------------------------------------------------
# Setup
#----------------------------------------------------------------
with open('../data/cleaned_data/cleaned_data.pkl', 'rb') as f:
	data = pickle.load(f)
feature_names = data['feature_names']

#----------------------------------------------------------------
# Methods and settings
#----------------------------------------------------------------
# method_list = ['pfs', 'pc_stable', 'si_hiton_pc', 'bnwsl', 'gfcsl', 'gfcl', 'glasso']
method_list = ['pfs', 'hpc_local', 'iamb_local', 'mmpc_local', 'pc_stable', 'gfcsl', 'glasso']
radius = 3
cluster_radii = [1, 2]
# p-value threshold (values greater than this will not be shown in table)
pval_threshold = 0.1

# silggm methods fdr
fdr = 0.05

# table[method][anchor][radius][system] = p_value
# node_counts[method][anchor][radius] = number of nodes in cluster (all nodes, not just cortical)
table = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
node_counts = defaultdict(lambda: defaultdict(dict))

#----------------------------------------------------------------
# DK → Yeo lookup (once)
#----------------------------------------------------------------
rows = []
with open('../data/dk_to_yeo.csv') as f:
	for row in csv.reader(f):
		if row and not row[0].startswith('#'):
			rows.append(row)

header = rows[0]
data_rows = rows[1:]
yeo_columns = header[1:header.index('TotalVertices')]

dk_to_yeo = {}
for row in data_rows:
	region = row[0]
	if region in ('TOTAL', 'NETWORK_p_value'):
		continue
	counts = {net: int(row[header.index(net)]) for net in yeo_columns}
	dk_to_yeo[region.replace('_', '').lower()] = max(counts, key=counts.get)

#----------------------------------------------------------------
# Main loop
#----------------------------------------------------------------
for method in method_list:

	if 'pfs' in method or method in ['hpc_local']:
		anchor_nodes = [f'Left\nCaudalmiddlefrontal\n(thick)']
	elif method in ['mmpc_local']:
		anchor_nodes = [f'Right\nPericalcarine\n(thick)']
	elif method in ['iamb_local', 'pc_stable', 'si_hiton_pc']:
		anchor_nodes = [f'Left\nPericalcarine\n(area)']
	elif method in ['bnwsl', 'gfcl', 'gfcsl']:
		anchor_nodes = [f'Left\nSuperiorparietal\n(thick)']
	elif method == 'glasso':
		anchor_nodes = [f'Left\nPericalcarine\n(area)']
	else:
		anchor_nodes = []

	anchor_node_names = [b.replace('\n', ' ') for b in anchor_nodes]

	with open(f'../results/hcp_{method}.pkl', 'rb') as f:
		result = pickle.load(f)

	target_features = result['target_features']

	if 'pfs' in method:
		A = prune_graph(
			result['Q'],
			target_features,
			result['qpath_max'],
			[0.15, 0.03, 0.02, 0.02],
			max_radius=radius
		)
	elif method in silggm_methods:
		A = result['adjacency_matrix'][fdr]
	else:
		A = result['adjacency_matrix']

	roi_to_system = build_roi_to_system(feature_names, dk_to_yeo)
	clean_feature_names = [format_feature_names(name) for name in feature_names]

	for cluster_radius in cluster_radii:
		print('################################')
		print(f'Radius: {cluster_radius+1}')
		print('################################')

		for i, anchor_node in enumerate(anchor_nodes):

			cluster = node_cluster(A, anchor_node=anchor_node, target_features=target_features,
				feature_names=clean_feature_names, max_radius=cluster_radius)

			node_counts[method][anchor_node_names[i]][cluster_radius] = len(cluster)

			component = []
			for f in (feature_names[j] for j in cluster):
				if not is_cortical_dk_feature(f):
					continue
				parts = f.split('_')
				hemi = 'Left' if parts[1] == 'L' else 'Right'
				component.append(f'{hemi} {parts[2].lower()}')

			enrich = yeo_enrichment(component, roi_to_system)

			print(f'{method} enrichment: {anchor_node_names[i]}')
			print('----------------------------------------------------------------')
			for system, stats in enrich.items():
				k = stats['overlap']
				K = stats['total']
				pval = stats['p_value']
				print(f'* {system}: {k}/{K} (p = {pval:.5f})')
			print()

			for system in yeo_columns:
				if system in enrich:
					table[method][anchor_node_names[i]][cluster_radius][system] = enrich[system]['p_value']
				else:
					table[method][anchor_node_names[i]][cluster_radius][system] = None

#----------------------------------------------------------------
# LaTeX table
#----------------------------------------------------------------
def latex_escape(s):
	return s.replace('_', r'\_')

system_names = {
	'Visual': 'VIS',
	'Somatomotor': 'SM',
	'DorsalAttention': 'DAN',
	'VentralAttention': 'VAN',
	'Frontoparietal': 'FPCN',
	'Default': 'DMN',
	'Limbic': 'LIM'
}

caption = (
    r'\caption{'
    r'\textit{Functional system enrichment of local brain networks associated with fluid cognition}. '
    r'For each method, we report $p$-values from hypergeometric over-representation analyses '
    r'for Yeo-7 functional networks~\cite{hcp_yeo} among cortical regions in the estimated local graph. '
    r'Local clusters are anchored at a single region of interest that connects directly to the target; '
    r'radii 2 and 3 denote cumulative neighborhoods containing all nodes within graph distances 2 '
    r'and 3 from the target that extend beyond this anchor. '
    r'The ``Nodes" column reports the total number of nodes in each cumulative cluster. '
    r'Columns correspond to Yeo-7 networks: VIS (visual), SM (somatomotor), DAN (dorsal attention), '
    r'VAN (ventral attention), LIM (limbic), FPCN (frontoparietal control), and DMN (default mode). '
    r'Only $p$-values less than $0.1$ are shown; entries marked ``--" indicate no enrichment at this threshold.'
    r'}'
)



systems_abbrev = [system_names[s] for s in yeo_columns]

def format_pval(p, threshold, eps=1e-4):
	if p is None or p > threshold:
		return '--'
	if p < eps:
		return r'$<10^{-4}$'
	return f'{p:.4f}'

col_spec = 'lcc' + 'c' * len(yeo_columns)
header = ['Method', 'Radius', 'Nodes'] + systems_abbrev

print(r'\begin{table*}[ht!]')
print(r'\centering')
print(r'\begin{tabular}{' + col_spec + '}')
print(r'\toprule')
print(' & '.join(header) + r' \\')
print(r'\bottomrule')

for method, anchors in table.items():
	for anchor, radii in anchors.items():

		if method not in method_names:
			method_names[method] = method

		method_label = latex_escape(method_names[method])
		radii_sorted = sorted(radii.keys())

		first_r = radii_sorted[0]
		row = [rf'\multirow{{{len(radii_sorted)}}}{{*}}{{{method_label}}}', str(first_r+1), str(node_counts[method][anchor][first_r])]
		for s in yeo_columns:
			row.append(format_pval(radii[first_r].get(s), pval_threshold))
		print(' & '.join(row) + r' \\')

		for r in radii_sorted[1:]:
			row = ['', str(r+1), str(node_counts[method][anchor][r])]
			for s in yeo_columns:
				row.append(format_pval(radii[r].get(s), pval_threshold))
			print(' & '.join(row) + r' \\')

		print(r'\hline')

print(r'\end{tabular}')
print(caption)
print(r'\label{tab:hcp}')
print(r'\end{table*}')
print()


