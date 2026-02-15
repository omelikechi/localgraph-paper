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

random_seed_list = [11201959]

#--------------------------------
# Load data
#--------------------------------
data_path = 'data/cleaned_data/cleaned_data.csv'
df = pd.read_csv(data_path)

target_names = ['histological_type', 'pathologic_stage', 'status']

# Identify target columns
target_features = [df.columns.get_loc(name) for name in target_names]

# Convert to numpy
X = df.to_numpy()
feature_names = df.columns.tolist()

#--------------------------------
# PFS arguments
#--------------------------------
radius = 1
qpath_max = 0.5
fdr_local = [0.5, 0.05, 0.0175]

# Gene-only neighborhoods
custom_nbhd = {}

#--------------------------------
# Prepare result dict
#--------------------------------
result = {
	'random_seed_list': random_seed_list,
	'qpath_max': qpath_max,
	'fdr_local': fdr_local,
	'custom_nbhd': custom_nbhd
}

# Add small noise to discrete targets
for idx in target_features:
	X[:, idx] += np.random.normal(0, 0.01, size=X.shape[0])

n, p = X.shape
result.update({
	'n': n,
	'p': p,
	'feature_names': feature_names,
	'target_features': target_features
})

#--------------------------------
# Run PFS
#--------------------------------
ipss_args = {'B': 200, 'preselector_args': {'expansion_factor': 2.5}}

for random_seed in random_seed_list:

	np.random.seed(random_seed)

	print(f'Random seed: {random_seed}')
	print(f'--------------------------------')

	result['A'] = pfs(
		X,
		target_features,
		qpath_max=qpath_max,
		fdr_local=fdr_local,
		custom_nbhd=custom_nbhd,
		feature_names=feature_names,
		max_radius=radius,
		ipss_args=ipss_args,
		verbose=True
	)

	#--------------------------------
	# Save result
	#--------------------------------
	if save_result:
		file_name = f'bc_pfs_r{radius}_{random_seed}'
		with open(f'{file_name}.pkl', 'wb') as f:
			pickle.dump(result, f)

	#--------------------------------
	# Plot result
	#--------------------------------
	if show_result:
		import matplotlib.pyplot as plt
		from localgraph import plot_graph

		rename_map = {
			'pathologic_stage': 'Stage',
			'status': 'Status',
			'histological_type': 'Histological\n type'
		}

		for old, new in rename_map.items():
			if old in feature_names:
				idx = feature_names.index(old)
				feature_names[idx] = new

		fig, ax = plt.subplots(1, 1, figsize=(17, 6))
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


