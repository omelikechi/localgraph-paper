# Plot results of Alzherimer's single cell data analyses

import os
import pickle

from localgraph import node_cluster, plot_graph, prune_graph
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import networkx as nx

import sys, os
sys.path.insert(0, os.path.abspath('../../..'))
from methods import method_names, silggm_methods

"""
Notes
----------------------------------------------------------------
- Glasso timed out after 24 hours (microglia)
- MMPC local timed out after 24 hours (microglia)
	- Took over 5 hours to find Markov blanket of certain nodes
- HPC local timed out after 24 hours (microglia)
	- Took several hours on multiple nodes
"""

#----------------------------------------------------------------
# Settings
#----------------------------------------------------------------
save_fig = True
dpi = 300

method = 'pfs'
cell_type = 'mg'

file_name = f'{method}_ad_{cell_type}'

# silggm methods fdr
fdr = 0.05

result_path = f'../results/{file_name}.pkl'

# layout and plot settings
radius = 3
node_size = 2000
font_size = 7
figsize = (18,9.5)
edge_digits = 2
show_weights = True if 'pfs' in method else False

# full cell type name
cell_type_name = {
	'astro':'Astrocytes',
	'mg':'Microglia',
	'OPC':'Oligodendrocyte progenitor cells'
}

#----------------------------------------------------------------
# Load result
#----------------------------------------------------------------
with open(result_path, 'rb') as f:
	result = pickle.load(f)

target_features = result['target_features']
feature_names = result['feature_names']

feature_names[-1] = 'AD status'

for key, value in result.items():
	if key not in ['random_seed', 'feature_names', 'target_features', 'Q', 'adjacency_matrix']:
		print(f'{key}: {value}')

if 'runtime' in result.keys():
	print(f'{method} runtime = {result['runtime']:.2f} seconds')
if 'A' in result.keys():
	result['adjacency_matrix'] = result['A']

# # format feature names
# for i, name in enumerate(feature_names):
# 	feature_names[i] = format_feature_names(name)

#----------------------------------------------------------------
# Load graph
#----------------------------------------------------------------
if method in silggm_methods:
	A = result['adjacency_matrix'][fdr]
elif 'pfs' in method:
	qpath_max = result['qpath_max']
	fdr_local = result['fdr_local']

	Q = result['Q']
	A = prune_graph(Q, target_features, qpath_max, fdr_local, max_radius=radius)

	# print radius 1 features
	print(f'\nCell type: {cell_type}')
	print(f'--------------------------------')
	for (i,j), q in Q.items():
		if i in target_features:
			 print(f'{feature_names[j]}: {q:.4f}')
else:
	A = result['adjacency_matrix']

# count nodes
nodes = set()
for i, j in A.keys():
	nodes.add(i)
	nodes.add(j)
n_nodes = len(nodes)
print(f'\nNumber of nodes in graph: {n_nodes}')

#----------------------------------------------------------------
# Plot
#----------------------------------------------------------------
plot_args = {
	'target_features': target_features,
	'feature_names': feature_names,
	'radius': radius,
	'include_outer_edges': False,
	'node_size': node_size,
	'font_size': font_size,
	'graph_layout': 'kk',
	'show_weights': show_weights,
	'edge_digits': edge_digits,
	'figsize': figsize,
	'dpi': dpi
}

fig, ax = plt.subplots(figsize=figsize)

plot_graph(graph=A, ax=ax, **plot_args)

# title
if method not in method_names:
	method_names[method] = method

ax.text(
	0.025, 0.965, f'{method_names[method]}: {cell_type_name[cell_type]}',
	transform=ax.transAxes,
	ha='left', va='top',
	fontsize=20, fontweight='bold'
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




