# Gene set over-representation analysis (snRNA-seq Alzheimer's data)
"""
Gene sets downloaded from https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H

- c2.cp.reactome.v2025.1.Hs.symbols.gmt (C2: curated gene sets)
	- Description: Gene sets in this collection are curated from various sources, including 
		online pathway databases and the biomedical literature. Reactome gene sets are derived 
		from Reactome and have been filtered to remove inter-set redundancy (see MSigDB release 
		notes for the current included Reactome version). http://www.reactome.org. For details 
		on citing gene sets in this collection see: http://www.reactome.org/cite.

- c5.go.bp.v2025.1.Hs.symbols.gmt (C5: ontology gene sets)
	- Description: Gene sets in this collection are derived from ontology resources. divided 
		into four sub collections derived from ontology annotations. Ontology annotations were 
		curated from databases maintained by their respective authorities. The C5:GO subcollection 
		is divided into three compoents (BP, CC, and MF) derived from Gene Ontology (GO), and 
		represent GO terms belonging to one of the three root GO ontologies: biological process 
		(BP), cellular component (CC), or molecular function (MF) respectively.
"""

import gseapy as gp
from localgraph import node_cluster, prune_graph
import pandas as pd
import pickle
import re

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 0)
pd.set_option('display.expand_frame_repr', False)

#----------------------------------------------------------------
# Import data and results
#----------------------------------------------------------------
method = 'pfs'

# p-value threshold and minimum overlap
p_thresh = 0.1
min_overlap = 5

cell_types = ['astro', 'mg', 'OPC']

for cell_type in cell_types:
	with open(f'../results/{method}_ad_{cell_type}.pkl', 'rb') as f:
		result = pickle.load(f)

	for key, item in result.items():
		if key not in ['feature_names', 'Q']:
			print(f'{key}: {item}')
	print()

for cell_type in cell_types:
	with open(f'../results/{method}_ad_{cell_type}.pkl', 'rb') as f:
		result = pickle.load(f)

	feature_names = result['feature_names']
	target_features = result['target_features']

	if method == 'pfs':
		Q = result['Q']
		radius = result['radius']
		qpath_max = result['qpath_max']
		fdr_local = result['fdr_local']
		Q = prune_graph(Q, target_features, qpath_max, fdr_local, max_radius=radius)
	else:
		Q = result['adjacency_matrix']

	#----------------------------------------------------------------
	# Genes
	#----------------------------------------------------------------
	anchor_node = 'ad_status'

	remove_targets = False if anchor_node == 'ad_status' else True
	cluster = node_cluster(Q, anchor_node=anchor_node, target_features=target_features, 
		feature_names=feature_names, max_radius=3, remove_targets=remove_targets)
	cluster = [feature_names[i] for i in cluster]
	if 'ad_status' in cluster:
		cluster.remove('ad_status')

	#----------------------------------------------------------------
	# Over-representation analysis
	#----------------------------------------------------------------
	# background set (all genes in analyzed dataset)
	background_set = feature_names

	selected_genes = sorted(set(cluster))

	gene_sets = ['c2.cp.reactome.v2025.1.Hs.symbols.gmt', 'c5.go.bp.v2025.1.Hs.symbols.gmt']

	for gene_set in gene_sets:
		file_name = f'../data/gene_sets/{gene_set}'

		print(f'\nGene set: {gene_set}')

		enr = gp.enrich(gene_list=selected_genes, gene_sets=file_name, background=background_set, outdir=None)

		# check if no results returned
		if not isinstance(enr.results, pd.DataFrame) or enr.results.empty:
			print('No enrichment results returned.')
			exit()

		res = enr.results[['Term', 'Overlap', 'Adjusted P-value', 'Genes']]
		# overlap and p-value filters
		overlap_ok = res['Overlap'].str.split('/').str[0].astype(int) >= min_overlap
		p_ok = res['Adjusted P-value'] < p_thresh
		res = res.loc[overlap_ok & p_ok].sort_values('Adjusted P-value')

		print_genes = False
		print(f'Cell type: {cell_type}')
		print(f'Size of cluster: {len(set(cluster))}')
		print(f'Size of background: {len(background_set)}')
		print(f'Minimum overlap: {min_overlap}')
		print(f'----------------------------------------------------------------')

		if print_genes:
			print(f'Pathway genes: {selected_genes}')

		def trunc(s,n):
			return s[:n] + '...' if isinstance(s, str) and len(s) > n else s

		res_to_print = res if print_genes else res.drop(columns=['Genes'])
		res_to_print = res_to_print.copy()
		res_to_print['Term'] = res_to_print['Term'].map(lambda x: trunc(x,50))

		print(res_to_print.head(10).to_string(index=False))
		print(f'----------------------------------------------------------------')


