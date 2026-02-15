# Apply graph estimation methods to Alzheimer's disease data

import pickle
import re

from localgraph import max_cor_response, plot_graph, restrict_to_local_graph
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

import sys, os
sys.path.insert(0, os.path.abspath('../..'))
from methods import run_method, silggm_methods

#----------------------------------------------------------------
# Setup
#----------------------------------------------------------------
save_result = False
show_result = True

random_seed = 6261928
np.random.seed(random_seed)

# cell types: astro, doublet, endo, mg, neuron, oligo, OPC, unID 
cell_type = 'OPC'

methods_to_run = ['mmpc_local']
apply_npn = False
max_cor_lambda = False
max_radius = 2
target_fdrs = [0.05]
bnlearn_fdr = 0.01

huge_methods = ['glasso', 'mb']

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
# Loop through methods
#----------------------------------------------------------------
for method_name in methods_to_run:

	print(f'Starting {method_name}')

	#----------------------------------------------------------------
	# bnlearn methods
	#----------------------------------------------------------------
	bnlearn_test = 'mi-g'
	bnlearn_methods = {
		'aracne': {'mi': 'mi-g'},
		'fast_iamb': {'alpha': bnlearn_fdr, 'test': bnlearn_test},
		'hpc': {'alpha': bnlearn_fdr, 'test': bnlearn_test},
		'iamb': {'alpha': bnlearn_fdr, 'test': bnlearn_test},
		'mmpc': {'alpha': bnlearn_fdr, 'test': bnlearn_test},
		'pc_stable': {'alpha': bnlearn_fdr, 'test': bnlearn_test},
		'si_hiton_pc': {'alpha': bnlearn_fdr, 'test': bnlearn_test}
	}

	if method_name in bnlearn_methods:
		method_args = bnlearn_methods[method_name]
		result = run_method(method_name, X, **method_args)
		A = result['adjacency_matrix']
		if max_radius is not None:
			result['adjacency_matrix'] = restrict_to_local_graph(
				A, target_features, max_radius
			)
		runtime = result['runtime']
		print(f'Runtime: {runtime:.2f} seconds')

	#----------------------------------------------------------------
	# bnlearn local methods
	#----------------------------------------------------------------
	bnlearn_test = 'mi-g'
	bnlearn_methods = {
		'fast_iamb_local': {'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features, 'verbose':True},
		'hpc_local': {'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features, 'verbose':True},
		'iamb_local': {'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features, 'verbose':True},
		'mmpc_local': {'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features, 'verbose':True},
		'pc_stable_local': {'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features, 'verbose':True},
		'si_hiton_pc_local': {'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features, 'verbose':True}
	}

	if method_name in bnlearn_methods:
		method_args = bnlearn_methods[method_name]
		result = run_method(method_name, X, **method_args)
		A = result['adjacency_matrix']
		if max_radius is not None:
			result['adjacency_matrix'] = restrict_to_local_graph(
				A, target_features, max_radius
			)
		runtime = result['runtime']
		print(f'Runtime: {runtime:.2f} seconds')

	#----------------------------------------------------------------
	# huge methods
	#----------------------------------------------------------------
	if method_name in huge_methods:
		X = StandardScaler().fit_transform(X)
		if max_cor_lambda:
			max_cors = max_cor_response(X, target_features)
			lambda_ = 0.99 * min(max_cors)
		else:
			lambda_ = None
		method_args = {'lambda_': lambda_, 'apply_npn': apply_npn}
		result = run_method(method_name, X, **method_args)
		A = result['adjacency_matrix']
		if max_radius is not None:
			result['adjacency_matrix'] = restrict_to_local_graph(
				A, target_features, max_radius
			)
		runtime = result['runtime']
		print(f'Runtime: {runtime:.2f} seconds')

	#----------------------------------------------------------------
	# Mixed graphical model
	#----------------------------------------------------------------
	if method_name == 'mgm':
		X = StandardScaler().fit_transform(X)
		method_args = {'cat_threshold': 4}
		result = run_method(method_name, X, **method_args)
		A = result['adjacency_matrix']
		if max_radius is not None:
			result['adjacency_matrix'] = restrict_to_local_graph(
				A, target_features, max_radius
			)
		runtime = result['runtime']
		print(f'Runtime: {runtime:.2f} seconds')

	#----------------------------------------------------------------
	# SILGGM
	#----------------------------------------------------------------
	if method_name in silggm_methods:
		X = StandardScaler().fit_transform(X)
		method_args = {'alpha': target_fdrs, 'apply_npn': apply_npn}
		result = run_method(method_name, X, **method_args)
		adjacency_matrix = result['adjacency_matrix']

		if max_radius is not None:
			if len(target_fdrs) == 1:
				alpha = target_fdrs[0]
				A = adjacency_matrix
				result['adjacency_matrix'] = {
					alpha: restrict_to_local_graph(A, target_features, max_radius)
				}
			else:
				for alpha in target_fdrs:
					A = adjacency_matrix[alpha]
					result['adjacency_matrix'][alpha] = restrict_to_local_graph(
						A, target_features, max_radius
					)
		runtime = result['runtime']
		print(f'Runtime: {runtime:.2f} seconds')

	#----------------------------------------------------------------
	# Plot result
	#----------------------------------------------------------------
	if show_result:
		plot_args = {
			'target_features': target_features,
			'feature_names': feature_names,
			'radius': max_radius,
			'font_size': 9,
			'node_size': 1800,
			'edge_widths': 3,
			'show_weights': False
		}

		if method_name not in silggm_methods:
			Q = result['adjacency_matrix']
			plot_graph(graph=Q, **plot_args)
		else:
			for target_fdr in target_fdrs:
				Q = result['adjacency_matrix'][target_fdr]
				print(f'Target FDR = {target_fdr}')
				plot_graph(graph=Q, **plot_args)

	#----------------------------------------------------------------
	# Save result
	#----------------------------------------------------------------
	if save_result:
		filename = f"ad_{method_name}_npn" if apply_npn else f"ad_{method_name}"
		if method_name in silggm_methods and len(target_fdrs) == 1:
			filename += f'_fdr{target_fdrs[0]}'
		filename += '.pkl'

		result['method_args'] = method_args
		result['target_features'] = target_features
		result['feature_names'] = feature_names
		result['apply_npn'] = apply_npn

		with open(filename, "wb") as f:
			pickle.dump(result, f)

		print(f"Saved result to {filename}")

	print(f'Finished {method_name}')



