# Simulation illustration for comparing truth, PFS, graphical lasso, and nodewise lasso (Figures 1a–1d)

from localgraph import pfs, plot_graph, subgraph_within_radius, tp_and_fp
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

import sys, os
sys.path.insert(0, os.path.abspath('..'))
from methods import run_method
from simulate_block import block_graph

#--------------------------------
# Settings
#--------------------------------
save_plot = False
dpi = 300
fig_name = f'illustration_dpi{dpi}'

# Show graph comparisons
show_graphs = True

# Simulation settings
random_seed = 168
np.random.seed(random_seed)
n = 200
radius = 3

# eignevalue scaling
lmin = 0.01
lmax = 10

# Graph block configuration
block_sizes = [1, 2, 97]
block_degree = [0, 0, 3]
connector_degree = [2, 3]
block_magnitude = [1, 1, 1]
connector_magnitude = [1, 1]
sigma = 0

# FDR thresholds
qpath_max = 0.4
fdr_local = [0.4, 0.4, 0.4]

# Method names
methods = ['truth', 'pfs', 'glasso', 'mb']
method_titles = {
	'truth': '(a) Truth',
	'pfs': '(b) Pathwise feature selection',
	'glasso': '(c) Graphical lasso',
	'mb': '(d) Nodewise lasso',
	
	# bnlearn
	'aracne':'(d) ARACNE',
	'fast_iamb':'FastIAMB',
	'hpc':'HPC',
	'iamb':'IAMB',
	'iamb_fdr':'IAMBFDR',
	'mmpc':'MMPC',
	'pc_stable':'PC',
	'si_hiton_pc':'SIHPC',
}

# Plotting options
figsize = (16,10)
cols = 2
plot_args = {'node_size': 1000, 'font_size': 14, 'edge_font_size': 9, 'edge_widths': 3}

#--------------------------------
# Main
#--------------------------------
# Generate data
data = block_graph(
	n=n, 
	lmin=lmin, 
	lmax=lmax, 
	block_sizes=block_sizes, 
	block_degree=block_degree, 
	block_magnitude=block_magnitude,
	connector_degree=connector_degree, 
	connector_magnitude=connector_magnitude, 
	random_seed=random_seed,
	sigma=sigma
)

X = data['X']
A_true = data['A']
target_features = data['target_features']
feature_names = data['feature_names']
max_cor_response = data['max_cor_response']
p = X.shape[1]

# Method configurations
lambda_ = 0.9999 * np.max(max_cor_response)

method_configs = {
	# bnlearn
	'aracne':{},
	'fast_iamb':{},
	'hpc':{},
	'iamb':{},
	'iamb_fdr':{},
	'mmpc':{},
	'pc_stable':{},
	'si_hiton_pc':{},

	# huge
	'glasso':{'lambda_': lambda_},
	'mb':{},

	# pfs
	'pfs': {
		'target_features':target_features,
		'selector':'l1',
		'qpath_max':qpath_max,
		'max_radius':len(fdr_local),
		'fdr_local':fdr_local,
		'criterion':'forward'
	}
}

# Convert PFS output to matrix
def dict_to_matrix(graph_dict, p):
	A = np.zeros((p, p))
	for (i,j), q in graph_dict.items():
		A[i,j] = q
		A[j,i] = q
	return A

def run_all_methods(method_name, X, **kwargs):
	if method_name == 'pfs':
		Q = pfs(X, **kwargs)
		return dict_to_matrix(Q,p)
	else:
		return run_method(method_name, X, **kwargs)['adjacency_matrix']

# Run all methods
results = {}
print(f'Method comparison')
print(f'--------------------------------')
print(f'Method name:')
print(f' - Full graph: true positives, false positives, FDP')
print(f' - Local graph: true positives, false positives, FDP')

for method_name in methods:
	if method_name == 'truth':
		results[method_name] = A_true
	else:
		A = run_all_methods(method_name, X, **method_configs[method_name])
		results[method_name] = A
		tp_full, fp_full = tp_and_fp(A, A_true, target_features)
		tp, fp = tp_and_fp(A, A_true, target_features, radius)
		print(f'{method_name}:')
		print(f' - Full graph: {tp_full}, {fp_full}, {fp_full / max(fp_full + tp_full, 1):.2f}')
		print(f' - Local graph: {tp}, {fp}, {fp / max(fp + tp, 1):.2f}')

#--------------------------------
# Plot estimated graphs
#--------------------------------
if show_graphs:
	A_true_local = subgraph_within_radius(A_true, target_features, radius)
	rows = (len(methods) + cols - 1) // cols
	fig, axes = plt.subplots(rows, cols, figsize=figsize)
	axes = axes.flatten()

	for i, method_name in enumerate(methods):
		A = results[method_name]
		plot_args['show_weights'] = (method_name == 'pfs')
		plot_graph(A, target_features, radius, true_graph=A_true_local, ax=axes[i], feature_names=feature_names, **plot_args)
		axes[i].text(
			0.025, 0.965, method_titles[method_name],
			transform=axes[i].transAxes,
			ha='left', va='top',
			fontsize=18, fontweight='bold'
		)

	for ax in axes[len(results):]:
		ax.set_visible(False)

	for ax in axes:
		if ax.get_visible():
			rect = patches.Rectangle((0, 0), 1, 1, transform=ax.transAxes, linewidth=2, edgecolor='gray', facecolor='none')
			ax.add_patch(rect)

	fig.tight_layout(pad=0.5)
	if save_plot:
		plt.savefig(f"{fig_name}.png", dpi=dpi)
	plt.show()
