# Plot results of HCP data analyses

import os
import pickle

from localgraph import node_cluster, plot_graph, prune_graph
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import networkx as nx

import sys, os
sys.path.insert(0, os.path.abspath('../../..'))
from applications.hcp.data.format_feature_names import format_feature_names
from methods import method_names, silggm_methods

"""
Notes
----------------------------------------------------------------
- Method with region of interest (ROI) in radius 1:
	- pfs, si_hiton_pc, pc_stable, bnwsl, gfcl, gfcsl, glasso
- Methods with no ROI in radius 1:
	- mb, dsgl, dsnwsl, aracne, mmpc
- Method that took too long to run:
	- fast_iamb, hpc
"""

#----------------------------------------------------------------
# Settings
#----------------------------------------------------------------
plot_adjusted_layout = False
save_fig = False
save_graph = False
dpi = 300

global_methods = []

# allow for single method (string) or list of methods
method = ['glasso', 'gfcsl', 'pc_stable', 'hpc_local', 'iamb_local', 'mmpc_local']
# method = 'iamb_local'

if isinstance(method, str):
	methods_to_plot = [method]
else:
	methods_to_plot = method
n_methods = len(methods_to_plot)

# silggm methods fdr
fdr = 0.05

# layout and plot settings
radius = 3
node_size = 3000 if n_methods == 1 else 750
font_size = 6 if n_methods == 1 else 4
figsize = (10,9) if n_methods == 1 else (18,9.5)
edge_digits = 2

#----------------------------------------------------------------
# Create figure/subplots
#----------------------------------------------------------------
if n_methods == 1:
	fig, axes = plt.subplots(figsize=figsize)
	axes = [axes]
else:
	ncols = min(3, n_methods)
	nrows = int((n_methods + ncols - 1) / ncols)
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
	axes = axes.flatten()

#----------------------------------------------------------------
# Loop over methods
#----------------------------------------------------------------
for idx, method in enumerate(methods_to_plot):

	file_name = f'hcp_{method}'
	result_path = f'../results/{file_name}.pkl'
	ax = axes[idx]

	#----------------------------
	# Load result
	#----------------------------
	with open(result_path, 'rb') as f:
		result = pickle.load(f)

	print(f'Method: {method}')
	print(f'-'*32)
	for key, item in result.items():
		if key in ['method_args']:
			print(f'{key}: {item}')
	print()

	target_features = result['target_features']
	feature_names = result['feature_names']

	if 'A' in result.keys():
		result['adjacency_matrix'] = result['A']

	# format feature names
	for i, name in enumerate(feature_names):
		feature_names[i] = format_feature_names(name)

	#----------------------------
	# Load graph
	#----------------------------
	if method in silggm_methods:
		A = result['adjacency_matrix'][fdr]
		show_weights = False
	elif 'pfs' in method:
		qpath_max = result['qpath_max']
		fdr_local = [0.15, 0.03, 0.02, 0.02]
		Q = result['Q']
		A = prune_graph(Q, target_features, qpath_max, fdr_local, max_radius=radius)
		show_weights = True
	else:
		A = result['adjacency_matrix']
		show_weights = False

	#----------------------------
	# Plot
	#----------------------------
	plot_args = {
		'target_features':target_features,
		'feature_names':feature_names,
		'radius':radius,
		'node_size':node_size,
		'font_size':font_size,
		'graph_layout':'kk',
		'show_weights':show_weights,
		'edge_digits':edge_digits,
		'figsize':figsize,
		'dpi':dpi,
		'save_graph':save_graph,
		'graph_name':file_name
	}

	plot_graph(graph=A, ax=ax, **plot_args)

	# title
	if method not in method_names:
		method_names[method] = method

	ax.text(
		0.025, 0.965, f'{method_names[method]}',
		transform=ax.transAxes,
		ha='left', va='top',
		fontsize=18, fontweight='bold'
	)

	# panel outline
	rect = patches.Rectangle(
		(0, 0), 1, 1,
		transform=ax.transAxes,
		linewidth=1.2,
		edgecolor='black',
		facecolor='none'
	)
	ax.add_patch(rect)

# remove unused axes
for j in range(idx+1, len(axes)):
	fig.delaxes(axes[j])

plt.tight_layout()

if save_fig:
	if n_methods == 1:
		plt.savefig(f'hcp_{methods_to_plot[0]}.png', dpi=dpi, bbox_inches='tight')
	else:
		plt.savefig('hcp_other_methods.png', dpi=dpi, bbox_inches='tight')

plt.show()


