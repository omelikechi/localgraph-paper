# Create the "limits of global error control" figure (Figure S1 in the Supplement)

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

save_plot = False

def create_graph(k1, k3):
	G = nx.Graph()

	# Add labeled nodes
	labeled_nodes = ['X_1', 'X_2', 'X_3']
	G.add_nodes_from(labeled_nodes)

	# Add edges between labeled nodes
	G.add_edge('X_1', 'X_2', color='black')
	G.add_edge('X_2', 'X_3', color='red')

	# Add unlabeled nodes connected to X_1
	unlabeled_nodes_X1 = [f'U1_{i}' for i in range(k1)]
	for node in unlabeled_nodes_X1:
		G.add_node(node)
		G.add_edge('X_1', node, color='black')

	# Add unlabeled nodes connected to X_3
	unlabeled_nodes_X3 = [f'U3_{i}' for i in range(k3)]
	for node in unlabeled_nodes_X3:
		G.add_node(node)
		G.add_edge('X_3', node, color='red')

	return G, labeled_nodes, unlabeled_nodes_X1, unlabeled_nodes_X3

def draw_graphs(G1, labeled_nodes1, unlabeled_nodes_X1_1, unlabeled_nodes_X3_1,
				G2, labeled_nodes2, unlabeled_nodes_X1_2, unlabeled_nodes_X3_2, save_plot=False):
	# Create figure with 2 side-by-side subplots
	fig, ax = plt.subplots(1, 2, figsize=(18, 5))

	def draw_single_graph(G, labeled_nodes, unlabeled_nodes_X1, unlabeled_nodes_X3, ax):
		# Extract edge colors
		edge_colors = [G[u][v]['color'] for u, v in G.edges]

		# Define fixed positions for labeled nodes (horizontally aligned)
		pos = {
			'X_1': (-1, 0),
			'X_2': (0, 0),
			'X_3': (1, 0),
		}

		# Arrange the unlabeled nodes **around** X_1 and X_3 in a structured way
		k1 = len(unlabeled_nodes_X1)
		k3 = len(unlabeled_nodes_X3)

		# Place nodes around X_1 in a semicircle (left side)
		angles_X1 = np.linspace(np.pi / 4, 7 * np.pi / 4, k1)  # Between 45 and 315 degrees
		radius = 0.75
		for i, node in enumerate(unlabeled_nodes_X1):
			pos[node] = (-1 + radius * np.cos(angles_X1[i]), radius * np.sin(angles_X1[i]))

		# Place nodes around X_3 in a semicircle (right side)
		angles_X3 = np.linspace(3 * np.pi / 4, -3 * np.pi / 4, k3)  # Between 135 and -135 degrees
		for i, node in enumerate(unlabeled_nodes_X3):
			pos[node] = (1 + radius * np.cos(angles_X3[i]), radius * np.sin(angles_X3[i]))

		# Define LaTeX-style labels
		latex_labels = {node: f'$X_{{{node[-1]}}}$' for node in labeled_nodes}

		# Draw nodes
		nx.draw_networkx_nodes(G, pos, nodelist=labeled_nodes, node_color='white', edgecolors='black', node_size=2000, ax=ax)
		nx.draw_networkx_nodes(G, pos, nodelist=unlabeled_nodes_X1, node_color='white', edgecolors='black', node_size=1000, ax=ax)
		nx.draw_networkx_nodes(G, pos, nodelist=unlabeled_nodes_X3, node_color='white', edgecolors='black', node_size=1000, ax=ax)

		# Draw edges
		nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=5, ax=ax)

		# Draw LaTeX labels only for labeled nodes
		nx.draw_networkx_labels(G, pos, labels=latex_labels, font_size=20, font_weight='bold', ax=ax)

		# Remove axis border
		ax.set_xticks([])
		ax.set_yticks([])
		ax.set_frame_on(False)

		# Draw a gray box around the subplot
		rect = patches.Rectangle((0, 0), 1, 1, transform=ax.transAxes, linewidth=2, edgecolor='gray', facecolor='none')
		ax.add_patch(rect)

	# Draw both graphs in separate subplots
	draw_single_graph(G1, labeled_nodes1, unlabeled_nodes_X1_1, unlabeled_nodes_X3_1, ax[0])
	draw_single_graph(G2, labeled_nodes2, unlabeled_nodes_X1_2, unlabeled_nodes_X3_2, ax[1])

	# Reduce white space
	fig.tight_layout(pad=0.5)

	# Save plot if needed
	if save_plot:
		plt.savefig('limits_error.png', dpi=300)
	
	plt.show()

# Example usage: Create two graphs with different k1 and k3 values
G1, labeled_nodes1, unlabeled_nodes_X1_1, unlabeled_nodes_X3_1 = create_graph(k1=8, k3=0)
G2, labeled_nodes2, unlabeled_nodes_X1_2, unlabeled_nodes_X3_2 = create_graph(k1=0, k3=8)

# Draw side-by-side graphs with gray boxes around each panel
draw_graphs(G1, labeled_nodes1, unlabeled_nodes_X1_1, unlabeled_nodes_X3_1,
			G2, labeled_nodes2, unlabeled_nodes_X1_2, unlabeled_nodes_X3_2, save_plot)






