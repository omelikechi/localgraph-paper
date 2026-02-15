# Plot results of PFS applied to TCGA breast cancer data (Figure 4 in the paper)

import pickle

from localgraph import plot_graph, prune_graph
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import networkx as nx

import sys, os
sys.path.insert(0, os.path.abspath('../../..'))
from applications.hcp.data.format_feature_names import format_feature_names

#--------------------------------
# Settings
#--------------------------------
save_fig = True
save_graph = False
dpi = 300

result_path = f'../results/hcp_pfs.pkl'
adjusted_path = 'hcp_pfs_adjusted.graphml'
graph_name = f'hcp_pfs'

# layout and plot settings
radius = 3
node_size = 3000
font_size = 6
figsize = (10,9)
edge_digits = 2

#--------------------------------
# Load result
#--------------------------------
with open(result_path, "rb") as f:
	result = pickle.load(f)

target_features = result['target_features']
feature_names = result['feature_names']

# format feature names
for i, name in enumerate(feature_names):
	feature_names[i] = format_feature_names(name)

#--------------------------------
# Load graph and layout
#--------------------------------
G = nx.read_graphml(adjusted_path)
G = nx.relabel_nodes(G, lambda x: int(x))
pos = {
	int(node): (float(data['x']), float(data['y']))
	for node, data in G.nodes(data=True)
}
A = None
show_weights = True

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
	'save_fig': save_fig,
	'graph_name': graph_name,
	'save_graph':save_graph,
	'dpi': dpi,
	'pos': pos
}

fig, ax = plt.subplots(figsize=figsize)

plot_graph(graph=G, ax=ax, **plot_args)

ax.text(
	0.025, 0.965, f'PFS',
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
	plt.savefig(f'{graph_name}.png', dpi=dpi, bbox_inches='tight')

plt.show()






