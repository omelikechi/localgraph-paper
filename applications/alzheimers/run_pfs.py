# Apply PFS to single cell Alzheimer's disease dataset
"""
Resources:
--------------------------------
- Data downloaded from https://adsn.ddnetbio.com
- Zenodo: https://zenodo.org/records/17302976
- Paper: "A single-cell atlas of entorhinal cortex from individuals with Alzheimer’s 
	disease reveals cell-type-specific gene expression regulation" by Grubman et al
- DOI: doi:10.1038/s41593-019-0539-4
""" 

import numpy as np
import pickle
import time
from localgraph import pfs

#----------------------------------------------------------------
# Setup
#----------------------------------------------------------------
save_result = True
show_result = True

random_seed = 6261928
np.random.seed(random_seed)

# cell types: astro, doublet, endo, mg, neuron, oligo, OPC, unID 
cell_type = 'astro'

# pfs parameters
radius = 1
qpath_max = 0.1
fdr_local = [0.025, 0.025, 0.03, 0.05]

#----------------------------------------------------------------
# Load data
#----------------------------------------------------------------
with open('./data/cleaned_data/cleaned_data.pkl', 'rb') as f:
	data = pickle.load(f)

X = data['X']
genes = data['genes']
meta = data['meta']

# cell type (ct)
mask = (meta['cellType'] == cell_type).values
X_ct = X[mask, :].toarray()
meta_ct = meta.loc[mask]
feature_names = list(genes)

# remove zero-variance genes
var = X_ct.var(axis=0)
keep_var = var > 1e-8
X_ct = X_ct[:, keep_var]
feature_names = [feature_names[i] for i in np.where(keep_var)[0]]

# add target
Y = (meta_ct['batchCond'] == 'AD').astype(int).values.reshape(-1,1)
X = np.hstack([X_ct, Y])
feature_names.append('ad_status')
target_features = [len(feature_names) - 1]

print(f'Cell type: {cell_type}')
print(f'Samples: {X.shape[0]}')
print(f'Genes: {X.shape[1]}\n')

#----------------------------------------------------------------
# Run PFS
#----------------------------------------------------------------
start = time.time()
Q = pfs(
	X,
	target_features,
	qpath_max=qpath_max,
	fdr_local=fdr_local,
	feature_names=feature_names,
	max_radius=radius,
	verbose=True
)

print('\nAD-associated local genes')
print('--------------------------------')
for (i,j), q in Q.items():
	if i == target_features[0] and q <= 0.25:
		print(f'{feature_names[j]}: {q:.3f}')

print(f'\nRuntime: {time.time() - start:.1f}')

#----------------------------------------------------------------
# Plot result
#----------------------------------------------------------------
if show_result:
	import matplotlib.pyplot as plt
	from localgraph import plot_graph

	fig, ax = plt.subplots(1, 1, figsize=(17,6))
	plot_graph(
		Q,
		ax=ax,
		target_features=target_features,
		feature_names=feature_names,
		radius=radius,
		node_size=1600,
		show_weights=True,
		font_size=8,
		graph_layout='kk',
		edge_digits=4
	)
	fig.tight_layout()
	plt.show()

#----------------------------------------------------------------
# Save results
#----------------------------------------------------------------
if save_result:

	result = {
		'cell_type': cell_type,
		'random_seed': random_seed,
		'radius': radius,
		'qpath_max': qpath_max,
		'fdr_local': fdr_local,
		'n': X_ct.shape[0],
		'p': X_ct.shape[1],
		'ad_cells': int(np.sum(Y)),
		'feature_names': feature_names,
		'target_features': target_features,
		'Q': Q
	}

	out_path = f'pfs_ad_{cell_type}.pkl'
	with open(out_path, 'wb') as f:
		pickle.dump(result, f)

	print(f'\nSaved result to {out_path}')


