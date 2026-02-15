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
save_graph = False
save_graph_settings = False
dpi = 300

result_path = f'../results/bc_pfs.pkl'
adjusted_path = '../results/bc_pfs_adjusted.graphml'
graph_name = f'breast_cancer_pfs_new_dpi{dpi}'

# layout and plot settings
radius = 3
node_size = 1000
font_size = 6
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

	with open('../results/bc_pfs_graph_settings.pkl', 'rb') as f:
		pfs_graph_settings = pickle.load(f)

	print(pfs_graph_settings['qpath_max'])
	print(pfs_graph_settings['fdr_local'])
	print(pfs_graph_settings['intermodal_fdr'])
	print(pfs_graph_settings['custom_nbhd']['CDH1\n(gene)'])

else:
	show_weights = True
	qpath_max = 1
	fdr_local = [0.25, 0.025, 0.0175]
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
			custom_nbhd[feature_name] = {'miR':intermodal_fdr, 'let':intermodal_fdr, 'protein':intermodal_fdr, 'gene':0.025}

	Q = result['Q']
	A = prune_graph(Q, target_features, qpath_max, fdr_local, max_radius,
					custom_nbhd=custom_nbhd, feature_names=feature_names)

	pos = None
	G = None

	# save graph settings
	if save_graph_settings:
		pfs_graph_settings = {
			'radius': radius,
			'qpath_max': qpath_max if not plot_adjusted_layout else None,
			'fdr_local': fdr_local if not plot_adjusted_layout else result['fdr_local'],
			'intermodal_fdr': intermodal_fdr if not plot_adjusted_layout else result['intermodal_fdr'],
			'custom_nbhd': custom_nbhd if not plot_adjusted_layout else None,
			'target_features': target_features,
			'feature_names': feature_names,
			'result_path': result_path,
		}

		settings_path = '../results/bc_pfs_graph_settings_new.pkl'
		with open(settings_path, 'wb') as f:
			pickle.dump(pfs_graph_settings, f)

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
	'save_graph':save_graph,
	'dpi': dpi,
	'pos': pos
}

# print features that appear in the plotted graph
graph_used = G if plot_adjusted_layout else A

plot_graph(graph=G if plot_adjusted_layout else A, **plot_args)
