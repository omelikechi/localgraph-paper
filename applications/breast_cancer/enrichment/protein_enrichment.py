# Compare estimated protein--protein edges to STRING database
"""
Source: https://string-db.org
"""

import pickle
import random

from localgraph import prune_graph
import numpy as np
import pandas as pd

#----------------------------------------------------------------
# Data
#----------------------------------------------------------------
data_path = '../data/cleaned_data/cleaned_data.csv'
df = pd.read_csv(data_path)

feature_names = df.columns.tolist()
idx_to_name = dict(enumerate(feature_names))

all_proteins = sorted({f.replace('#','') for f in feature_names if f.endswith('#')})

#----------------------------------------------------------------
# Methods
#----------------------------------------------------------------
methods_to_run = {
	'pfs':'PFS',
	'aracne':'ARACNE',
	'gfcsl':'GFCSL',
	'gfcsl_npn':'GFCSL(NPN)',
	'dsnwsl':'DSNWSL',
	'dsnwsl_npn':'DSNWSL(NPN)',
	'mb':'NLasso',
	'mb_npn':'NLasso(NPN)'
}

B = 100000
string_score_threshold = 400

#----------------------------------------------------------------
# load STRING (local)
#----------------------------------------------------------------
links_path = '../data/stringdb/9606.protein.links.v12.0.txt'
info_path = '../data/stringdb/9606.protein.info.v12.0.txt'

links = pd.read_csv(links_path, sep=' ')
info = pd.read_csv(info_path, sep='\t')[['#string_protein_id', 'preferred_name']]
info = info.rename(columns={'#string_protein_id': 'protein_id'})

links = links.merge(
	info, left_on='protein1', right_on='protein_id', how='left'
).rename(columns={'preferred_name': 'gene1'}).drop(columns='protein_id')

links = links.merge(
	info, left_on='protein2', right_on='protein_id', how='left'
).rename(columns={'preferred_name': 'gene2'}).drop(columns='protein_id')

links = links[['gene1', 'gene2', 'combined_score']].dropna()

string_dict = {
	tuple(sorted((row.gene1, row.gene2))): row.combined_score
	for row in links.itertuples(index=False)
}

def string_score(p1, p2):
	return string_dict.get(tuple(sorted((p1, p2))), 0)

#----------------------------------------------------------------
# Run all methods
#----------------------------------------------------------------
rows = []

for method in methods_to_run:

	# silggm-style FDR handling
	if '_npn' in method:
		fdr = 0.15 if method == 'gfcsl_npn' else 0.5
	else:
		fdr = 0.1 if method == 'gfcsl' else 0.5

	result_path = f'../results/bc_{method}.pkl'
	with open(result_path, 'rb') as f:
		result = pickle.load(f)

	if method == 'pfs':
		Q = result['Q']
		target_features = result['target_features']

		with open('../results/bc_pfs_graph_settings.pkl', 'rb') as f:
			pfs_graph_settings = pickle.load(f)

		A = prune_graph(
			Q,
			target_features,
			pfs_graph_settings['qpath_max'],
			pfs_graph_settings['fdr_local'],
			pfs_graph_settings['radius'],
			custom_nbhd=pfs_graph_settings['custom_nbhd'],
			feature_names=feature_names
		)

	elif method in ['dsnwsl', 'dsnwsl_npn', 'gfcsl', 'gfcsl_npn']:
		A = result['adjacency_matrix'][fdr]

	else:
		A = result['adjacency_matrix']

	inferred_edges = sorted({
		tuple(sorted((
			idx_to_name[i].replace('#', ''),
			idx_to_name[j].replace('#', '')
		)))
		for (i,j) in A
		if '#' in idx_to_name[i] and '#' in idx_to_name[j]
	})

	# total edges
	total_edges = len(inferred_edges)

	# observed
	observed_scores = np.array([string_score(p1, p2) for p1, p2 in inferred_edges])
	observed_supported = np.sum(observed_scores >= string_score_threshold)

	# null
	null_supported = []
	for _ in range(B):
		rand_edges = [tuple(random.sample(all_proteins, 2)) for _ in inferred_edges]
		scores = np.array([string_score(p1, p2) for p1, p2 in rand_edges])
		null_supported.append(np.sum(scores >= string_score_threshold))

	null_supported = np.array(null_supported)

	p_empirical = (np.sum(null_supported >= observed_supported) + 1) / (B + 1)

	rows.append({
		'Method': methods_to_run[method],
		'Total': total_edges,
		'Supported': observed_supported,
		'Expected': null_supported.mean(),
		'Fold enrichment': observed_supported / max(null_supported.mean(), 1e-6),
		r'Empirical $p$-value': p_empirical
	})


#----------------------------------------------------------------
# Final table
#----------------------------------------------------------------
results = pd.DataFrame(rows)[[
	'Method',
	'Total',
	'Supported',
	'Expected',
	'Fold enrichment',
	r'Empirical $p$-value'
]]

print(results)

print(
	results.to_latex(
		index=False,
		float_format='%.5g',
		column_format='lccccc'
	)
)


