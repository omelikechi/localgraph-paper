# Apply PFS to TCGA breast cancer data
"""
- Targets: histological type, pathologic stage, and status (alive or deceased)
- Covariates: genes (RNAseq), proteins (RPPA), and miRNAs
- Note: PFS is a stochastic algorithm, so results may vary slightly across machines
	even when using the same random seed.
"""

import pickle
import numpy as np
import pandas as pd

from localgraph import pfs

#--------------------------------
# Setup
#--------------------------------
save_result = False
show_result = True

random_seed = 11201959
np.random.seed(random_seed)

#--------------------------------
# Load data
#--------------------------------
data_path = 'data/cleaned_data/cleaned_data.csv'
df = pd.read_csv(data_path)

# Identify target columns
target_names = ['histological_type', 'pathologic_stage', 'status']
target_features = [df.columns.get_loc(name) for name in target_names]

# Convert data to numpy array
X = df.to_numpy()
feature_names = df.columns.tolist()

#--------------------------------
# PFS arguments
#--------------------------------
radius = 3
qpath_max = 1
fdr_local = [0.25, 0.023, 0.0175]

# Custom neighborhoods
custom_nbhd = {name:{'nbhd_fdr':0.2 if name != 'histological_type' else 0.25} for name in target_names}
intermodal_fdr = 0.035

for feature_name in feature_names:
	if feature_name in target_names:
		continue
	elif '#' in feature_name:
		custom_nbhd[feature_name] = {'miR':intermodal_fdr, 'let':intermodal_fdr, '*':intermodal_fdr}
	elif 'miR-' in feature_name or 'let-' in feature_name:
		custom_nbhd[feature_name] = {'#':intermodal_fdr, '*':intermodal_fdr}
	elif '*' in feature_name:
		custom_nbhd[feature_name] = {'miR':intermodal_fdr, 'let':intermodal_fdr, '#':intermodal_fdr}

#--------------------------------
# Prepare result dict
#--------------------------------
result = {
	'random_seed': random_seed,
	'qpath_max': qpath_max,
	'fdr_local': fdr_local,
	'custom_nbhd': custom_nbhd,
	'intermodal_fdr': intermodal_fdr
}

# Add noise to discrete target variables to reduce bias in tree-based importance scores
for idx in target_features:
	X[:, idx] += np.random.normal(0, 0.01, size=X.shape[0])

n, p = X.shape
result.update({'n': n, 'p': p, 'feature_names': feature_names, 'target_features': target_features})

#--------------------------------
# Run PFS
#--------------------------------
ipss_args = {'B':200, 'preselector_args':{'expansion_factor':2.5}}
result['A'] = pfs(X, target_features, qpath_max=qpath_max, fdr_local=fdr_local,
	custom_nbhd=custom_nbhd, feature_names=feature_names, max_radius=radius,
	ipss_args=ipss_args, verbose=True)

#--------------------------------
# Save result
#--------------------------------
if save_result:
	with open('breast_cancer_pfs_results.pkl', 'wb') as f:
		pickle.dump(result, f)

#--------------------------------
# Plot result
#--------------------------------
if show_result:
	import matplotlib.pyplot as plt
	from localgraph import plot_graph

	# Rename target nodes
	rename_map = {
		'pathologic_stage': 'Stage',
		'status': 'Status',
		'histological_type': 'Histological\n type'
	}

	for old, new in rename_map.items():
		if old in feature_names:
			idx = feature_names.index(old)
			feature_names[idx] = new

	# Plot
	fig, ax = plt.subplots(1, 1, figsize=(17,6))
	plot_graph(
		result['A'],
		ax=ax,
		target_features=target_features,
		feature_names=feature_names,
		radius=radius,
		node_size=1750,
		show_weights=True,
		font_size=9,
		graph_layout='kk',
		edge_digits=4
	)
	fig.tight_layout()
	plt.show()


