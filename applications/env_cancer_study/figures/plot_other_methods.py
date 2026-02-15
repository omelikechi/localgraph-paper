# Plot results of the graphical lasso applied to county-level cancer data (Figure 1e in the paper)

import os
import pickle

from localgraph import plot_graph, prune_graph
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import networkx as nx

#--------------------------------
# Settings
#--------------------------------
plot_adjusted_layout = False
save_fig = True
dpi = 300

method = 'pc_stable'
include_npn = False

file_name = f'eqi_{method}'

result_path = f'../results/{file_name}.pkl'

# layout and plot settings
radius = 2
node_size = 2000 if include_npn else 3500
font_size = 6 if include_npn else 8
figsize = (18,6) if include_npn else (16,9)
edge_digits = 3
show_weights = False

# method names in plot
method_names = {
	# bnlearn
	'aracne':'ARACNE',
	'hpc':'HPC',
	'mmpc':'MMPC',
	'pc_stable':'PC',
	'sihpc':'SIHPC',

	# huge
	'glasso':'Graphical lasso',
	'mb':'Nodewise lasso',

	# silggm
	'bnwsl':'BNWSL',
	'dsgl':'DSGL',
	'dsnwsl':'DSNWSL',
	'gfcl':'GFCL',
	'gfcsl':'GFCSL'
}

#--------------------------------
# Load result
#--------------------------------
with open(result_path, "rb") as f:
	result = pickle.load(f)

target_features = result['target_features']
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
# feature names
manual_feature_name_map = {
	'VacantUnits':'Vacant\nUnits',
	'cat_RUCC':'Rural',
	'Over21':'Age',
	'pm2.5(A)':f'PM$_{{2.5}}$(A)',
	'pm10(A)':f'PM$_{{10}}$(A)',
	'HomeValue':'Home\nValue',
	'Rooms/Home':'Rooms',
	'NaPct(L)':'Na(L)',
	'Unemployed':'Employed',
	'C3H3N(A)':'Acrylonitrile(A)', 
	'N2H2(A)':'Hydrazine(A)',
	'EOx(A)':'EtO(A)',
	'C2HCl3(A)':'TCE(A)',
	'ViolentCrime':'Violent\ncrime',
	'EntBusiness':'Leisure\nbusinesses',
	'fatal_rate_log':'Traffic\nfatalities',
	'WorkOutCo':'Employed out\nof county'
}

feature_names = [manual_feature_name_map.get(name, name) for name in feature_names]

# Apply suffix replacements for unmapped names
def replace_suffix(name):
	if name.endswith('(A)'):
		return name[:-3] + '\n(air)'
	elif name.endswith('(W)'):
		return name[:-3] + '\n(water)'
	elif name.endswith('(L)'):
		return name[:-3] + '\n(land)'
	else:
		return name
feature_names = [replace_suffix(name) for name in feature_names]

#--------------------------------
# Load graph and layout
#--------------------------------
max_radius = radius
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

	for ax in axes:
		for spine in ax.spines.values():
			spine.set_visible(True)
			spine.set_linewidth(0.8)
			spine.set_color('black')

	# without nonparanormal
	plot_graph(graph=result['adjacency_matrix'], ax=axes[0], **plot_args)

	# with nonparanormal
	plot_graph(graph=result_npn['adjacency_matrix'], ax=axes[1], **plot_args)
	
	# titles
	axes[0].text(0.025, 0.965, f'(a) {method_names[method]}', transform=axes[0].transAxes,
		ha='left', va='top', fontsize=14, fontweight='bold')
	axes[1].text(0.025, 0.965, f'(b) {method_names[method]} (NPN)', transform=axes[1].transAxes,
		ha='left', va='top', fontsize=14, fontweight='bold')

	# panel outlines
	for ax in axes:
		rect = patches.Rectangle((0, 0), 1, 1, transform=ax.transAxes, linewidth=1.2, edgecolor='black', facecolor='none')
		ax.add_patch(rect)

else:
	fig, ax = plt.subplots(figsize=figsize)
	plot_graph(graph=result['adjacency_matrix'], ax=ax, **plot_args)

	# title
	ax.text(0.025, 0.965, f'{method_names[method]}', transform=ax.transAxes,
		ha='left', va='top', fontsize=18, fontweight='bold')

	# panel outline
	rect = patches.Rectangle((0, 0), 1, 1, transform=ax.transAxes, linewidth=1.2, edgecolor='black', facecolor='none')
	ax.add_patch(rect)

plt.tight_layout()

if save_fig:
	plt.savefig(f'{file_name}.png', dpi=dpi, bbox_inches='tight')

plt.show()



