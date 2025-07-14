# Plot results of PFS applied to TCGA breast cancer data (Figure 4 in the paper)

import os
import pickle

from localgraph import plot_graph, prune_graph
import matplotlib.pyplot as plt
import networkx as nx

#--------------------------------
# Settings
#--------------------------------
plot_adjusted_layout = True
save_fig = False
dpi = 300

result_path = f'../results/breast_cancer_pfs_results.pkl'
adjusted_path = '../results/adjusted_graph.graphml'
graph_name = f'breast_cancer_pfs_dpi{dpi}'

# layout and plot settings
radius = 3
node_size = 1500
font_size = 7
figsize = (16,9)
edge_digits = 2

#--------------------------------
# Load result
#--------------------------------
with open(result_path, "rb") as f:
	result = pickle.load(f)

target_features = result['target_features']
feature_names = result['feature_names']

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
# Load graph and layout
#--------------------------------
if plot_adjusted_layout:
	G = nx.read_graphml(adjusted_path)
	G = nx.relabel_nodes(G, lambda x: int(x))
	pos = {
		int(node): (float(data['x']), float(data['y']))
		for node, data in G.nodes(data=True)
	}
	A = None
	show_weights = True
else:
	show_weights = True
	qpath_max = 1
	fdr_local = [0.3, 0.023, 0.0175]
	max_radius = radius
	custom_nbhd = {}
	intermodal_fdr = 0.035
	for feature_name in feature_names:
		if feature_name in custom_nbhd:
			continue
		elif 'protein' in feature_name:
			custom_nbhd[feature_name] = {'miR':intermodal_fdr, 'let':intermodal_fdr, 'gene':intermodal_fdr}
		elif 'miR-' in feature_name or 'let-' in feature_name:
			custom_nbhd[feature_name] = {'protein':intermodal_fdr, 'gene':intermodal_fdr}
		elif 'gene' in feature_name:
			custom_nbhd[feature_name] = {'miR':intermodal_fdr, 'let':intermodal_fdr, 'protein':intermodal_fdr}

	Q = result['Q']
	A = prune_graph(Q, target_features, qpath_max, fdr_local, max_radius,
					custom_nbhd=custom_nbhd, feature_names=feature_names)

	pos = None
	G = None

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
	'save_fig': save_fig,
	'graph_name': graph_name,
	'dpi': dpi,
	'pos': pos
}

plot_graph(graph=G if plot_adjusted_layout else A, **plot_args)
