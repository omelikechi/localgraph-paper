# Apply PFS to HCP neuroimaging data
"""
- Target: NIH Toolbox Fluid Cognition composite
- Covariates: structural morphometry + reduced, summary-level phenotypes
- Input: cleaned HCP Young Adult dataset (see data preparation script)
- Purpose: estimate local dependency structure around cognition using PFS
"""

import time

from collections import defaultdict
import pickle
import numpy as np

from localgraph import pfs

#--------------------------------
# Setup
#--------------------------------
save_result = False
show_result = True

random_seed = 1
np.random.seed(random_seed)

#----------------------------------------------------------------
# Load data
#----------------------------------------------------------------
with open('./data/cleaned_data/cleaned_data.pkl', 'rb') as f:
	data = pickle.load(f)

X = data['X']
feature_names = data['feature_names']
target_features = data['target_indices']
age_adjusted = data['age_adjusted']

n, p = X.shape

print(f'n = {n}')
print(f'p = {p}')

#----------------------------------------------------------------
# PFS arguments
#----------------------------------------------------------------
radius = 4
qpath_max = 0.25
fdr_local = [0.15, 0.03, 0.03, 0.05]

ipss_args = {'B': 200}

#----------------------------------------------------------------
# Prepare result container
#----------------------------------------------------------------
result = {
	'random_seed': random_seed,
	'age_adjusted': age_adjusted,
	'qpath_max': qpath_max,
	'fdr_local': fdr_local,
	'ipss_args': ipss_args,
	'n': n,
	'p': p,
	'feature_names': feature_names,
	'target_features': target_features
}

#----------------------------------------------------------------
# Run PFS
#----------------------------------------------------------------
start = time.time()
Q = pfs(
	X,
	target_features,
	qpath_max=qpath_max,
	fdr_local=fdr_local,
	feature_names=feature_names,
	max_radius=radius,
	ipss_args=ipss_args,
	verbose=True
)
runtime = time.time() - start

result['runtime'] = runtime
result['Q'] = Q

#----------------------------------------------------------------
# Print results by radius
#----------------------------------------------------------------
def print_nodes_by_radius(Q, target_features, feature_names, max_radius):

	adj = defaultdict(set)
	for (i, j) in Q.keys():
		adj[i].add(j)
		adj[j].add(i)

	for t in target_features:
		print(f'\nTarget: {feature_names[t]}')
		print('-' * 40)

		visited = {t}
		frontier = {t}

		for r in range(1, max_radius + 1):
			next_frontier = set()
			for u in frontier:
				next_frontier |= adj[u]

			next_frontier -= visited
			visited |= next_frontier

			if not next_frontier:
				break

			print(f'Radius {r}:')
			for v in sorted(next_frontier):
				q = min(
					Q.get((u, v), np.inf)
					for u in frontier
					if (u, v) in Q or (v, u) in Q
				)
				print(f'  {feature_names[v]} (q = {q:.3f})')

			frontier = next_frontier

print_nodes_by_radius(
	Q,
	target_features=target_features,
	feature_names=feature_names,
	max_radius=radius
)
print()

#----------------------------------------------------------------
# Plot result
#----------------------------------------------------------------
if show_result:
	import matplotlib.pyplot as plt
	from localgraph import plot_graph

	fig, ax = plt.subplots(1, 1, figsize=(17, 6))
	plot_graph(
		result['Q'],
		ax=ax,
		target_features=target_features,
		feature_names=feature_names,
		radius=radius,
		node_size=1600,
		show_weights=True,
		font_size=8,
		graph_layout='kk',
		edge_digits=4
	)
	fig.tight_layout()
	plt.show()

#----------------------------------------------------------------
# Save results
#----------------------------------------------------------------
if save_result:
	out_path = 'pfs_hcp_result_new.pkl'
	with open(out_path, 'wb') as f:
		pickle.dump(result, f)

	print(f'Saved PFS result to {out_path}')
