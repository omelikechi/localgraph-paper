# Plot results of PFS applied to county-level cancer data (Figure 3 in the paper)

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

result_path = f'../results/env_study_pfs_results.pkl'
adjusted_path = '../results/adjusted_graph_pfs.graphml'
graph_name = f'env_study_pfs_dpi{dpi}'

# layout and plot settings
radius = 2
node_size = 3000
font_size = 8
figsize = (16,9)
edge_digits = 3
show_weights = True

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
	'EntBusiness':'Leisure\nbusinesses'
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
if plot_adjusted_layout:
	G = nx.read_graphml(adjusted_path)
	G = nx.relabel_nodes(G, lambda x: int(x))  # convert node labels to int
	pos = {
		int(node): (float(data['x']), float(data['y']))
		for node, data in G.nodes(data=True)
	}
	A = None
	
else:
	qpath_max = 1
	fdr_local = [0.04, 0.008]
	max_radius = radius
	
	# custom neighborhood
	"""
	- Education: Many variables with the same q-value. Used efp scores to break ties. Custom variables ordered
		in terms of efp scores (Poverty: 0.04, Income: 0.04, HomeValue: 0.05, EntBusiness: 0.05, Colonoscopy: 0.06
		ViolentCrime: 0.06, White: 0.06, Mortality: 0.07)
	"""
	custom_nbhd = {
		'Mortality':{'nbhd_fdr':0.009},
		'Education':{'nbhd_fdr':0.001, 'Poverty':0.1, 'Income':0.1, 'HomeValue':0.1, 'EntBusiness':0.1, 
			'Colonoscopy':0.1, 'ViolentCrime':0.1, 'White':0.1, 'Mortality':0.1},
		'Hispanic':{'nbhd_fdr':0.01},
		'pm2.5(A)':{'Education':0.001},
		'Poverty':{'nbhd_fdr':0.011},
		'Smoking':{'nbhd_fdr':0.01}
	}

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
