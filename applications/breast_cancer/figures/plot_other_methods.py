# Plot results of other methods applied to TCGA breast cancer data

import os
import pickle

from localgraph import plot_graph, prune_graph
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import networkx as nx

# Notes
"""
- Graphical lasso did not complete within 24 hours
- Graphical lasso with nonparanormal did not complete withing 24 hours
"""

#--------------------------------
# Settings
#--------------------------------
plot_adjusted_layout = False
save_fig = False
dpi = 300

method = 'dsnwsl'

file_name = f'bc_{method}'

fdr = 0.1 if method == 'gfcsl' else 0.5
include_npn = True
npn_fdr = 0.15 if method == 'gfcsl' else 0.5

fdr_methods = ['bnwsl', 'dsgl', 'dsnwsl', 'gfcl', 'gfcsl']

result_path = f'../results/{file_name}.pkl'

# layout and plot settings
radius = 3
node_size = 800 if include_npn else 2000
font_size = 3 if include_npn else 8
figsize = (18, 6) if include_npn else (16,9)
edge_digits = 4
show_weights = True if method == 'pfs' else False

# method names in plot
method_names = {
	# bnlearn
	'aracne':'ARACNE',
	'hpc':'HPC',
	'mmpc':'MMPC',
	'pc_stable':'PC',
	'pfs':'PFS',
	'sihpc':'SIHPC',

	# huge
	'glasso':'Graphical lasso',
	'mb':'Nodewise lasso',

	# silggm
	'dsnwsl':'DSNWSL',
	'gfcsl':'GFCSL'
}

#--------------------------------
# Load result
#--------------------------------
with open(result_path, 'rb') as f:
	result = pickle.load(f)

target_features = result['target_features']
# target_features = [1]
feature_names = result['feature_names']

if 'runtime' in result.keys():
	print(f'{method} runtime = {result['runtime']:.2f} seconds')
if 'A' in result.keys():
	result['adjacency_matrix'] = result['A']

if include_npn:
	with open(f'../results/{file_name}_npn.pkl', 'rb') as f:
		result_npn = pickle.load(f)
	if 'A' in result_npn.keys():
		result_npn['adjacency_matrix'] = result_npn['A']

#--------------------------------
# Format feature names
#--------------------------------
rename_dict = {
	'histological_type': 'Histologic\ntype',
	'pathologic_stage': 'Stage',
	'status': 'Status'
}

for i, name in enumerate(feature_names):
	if name in rename_dict:
		feature_names[i] = rename_dict[name]
	elif '#' in name:
		feature_names[i] = name.replace('#', '\n(protein)')
	elif '*' in name:
		feature_names[i] = name.replace('*', '\n(gene)')

#--------------------------------
# Load graph
#--------------------------------
if method in fdr_methods:
	A = result['adjacency_matrix'][fdr]
else:
	A = result['adjacency_matrix']

#--------------------------------
# Plot
#--------------------------------
plot_args = {
	'target_features': target_features,
	'feature_names': feature_names,
	'radius': radius,
	'node_size': node_size,
	'font_size': font_size,
	'graph_layout': 'kk',
	'show_weights': show_weights,
	'edge_digits': edge_digits,
	'figsize': figsize,
	'dpi': dpi
}

if include_npn:
	fig, axes = plt.subplots(1, 2, figsize=figsize)

	# without NPN
	plot_graph(graph=result['adjacency_matrix'][fdr], ax=axes[0], **plot_args)

	# with NPN
	plot_graph(graph=result_npn['adjacency_matrix'][npn_fdr], ax=axes[1], **plot_args)

	# titles
	axes[0].text(
		0.025, 0.965, f'(a) {method_names[method]}',
		transform=axes[0].transAxes,
		ha='left', va='top',
		fontsize=14, fontweight='bold'
	)
	axes[1].text(
		0.025, 0.965, f'(b) {method_names[method]} (NPN)',
		transform=axes[1].transAxes,
		ha='left', va='top',
		fontsize=14, fontweight='bold'
	)

	# panel outlines
	for ax in axes:
		rect = patches.Rectangle(
			(0, 0), 1, 1,
			transform=ax.transAxes,
			linewidth=1.2,
			edgecolor='black',
			facecolor='none'
		)
		ax.add_patch(rect)

else:
	fig, ax = plt.subplots(figsize=figsize)

	plot_graph(graph=A, ax=ax, **plot_args)

	# title
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

plt.tight_layout()

if save_fig:
	plt.savefig(f'{file_name}.png', dpi=dpi, bbox_inches='tight')

plt.show()




