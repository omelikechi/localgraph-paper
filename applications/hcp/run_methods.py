# Apply graph estimation methods to HCP neuroimaging data
"""
- Targets: focal cortical thickness measures
- Covariates: structural morphometry + selected phenotypic variables
"""

import pickle
import re

from localgraph import max_cor_response, plot_graph, restrict_to_local_graph
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

import sys, os
sys.path.insert(0, os.path.abspath('../..'))
from methods import run_method, method_type, bnlearn_methods, bnlearn_local_methods, huge_methods, silggm_methods

#----------------------------------------------------------------
# Setup
#----------------------------------------------------------------
save_result = True
plot_result = True

random_seed = 302
np.random.seed(random_seed)

methods_to_run = ['mmpc_local']
apply_npn = False
max_cor_lambda = False
max_radius = 3
target_fdrs = [0.05]
bnlearn_fdr = 0.05
bnlearn_local_fdr = 0.025

#----------------------------------------------------------------
# Load data
#----------------------------------------------------------------
with open('./data/cleaned_data/cleaned_data.pkl', 'rb') as f:
	data = pickle.load(f)

X = data['X']
feature_names = data['feature_names']
target_features = data['target_indices']
age_adjusted = data['age_adjusted']
n = data['n']
p = data['p']

#----------------------------------------------------------------
# Loop through methods
#----------------------------------------------------------------
for method_name in methods_to_run:

	print(f'Starting {method_name}')

	mtype = method_type(method_name)

	print(mtype)

	#----------------------------------------------------------------
	# bnlearn methods
	#----------------------------------------------------------------
	if mtype == 'bnlearn':
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

		method_args = bnlearn_methods[method_name]
		result = run_method(method_name, X, **method_args)
		A = result['adjacency_matrix']
		if max_radius is not None:
			result['adjacency_matrix'] = restrict_to_local_graph(A, target_features, max_radius)
		runtime = result['runtime']
		print(f'Runtime: {runtime:.2f} seconds')

	#----------------------------------------------------------------
	# bnlearn local methods
	#----------------------------------------------------------------
	if mtype == 'bnlearn_local':
		bnlearn_test = 'mi-g'
		verbose = True
		bnlearn_local_methods = {
			'fast_iamb_local': {'alpha':bnlearn_local_fdr, 'test':bnlearn_test, 'target_features':target_features, 
				'radius':max_radius, 'verbose':verbose},
			'hpc_local': {'alpha':bnlearn_local_fdr, 'test':bnlearn_test, 'target_features':target_features, 
				'radius':max_radius, 'verbose':verbose},
			'iamb_local': {'alpha':bnlearn_local_fdr, 'test':bnlearn_test, 'target_features':target_features, 
				'radius':max_radius, 'verbose':verbose},
			'mmpc_local': {'alpha':bnlearn_local_fdr, 'test':bnlearn_test, 'target_features':target_features, 
				'radius':max_radius, 'verbose':verbose},
			'pc_stable_local': {'alpha':bnlearn_local_fdr, 'test':bnlearn_test, 'target_features':target_features, 
				'radius':max_radius, 'verbose':verbose},
			'si_hiton_pc_local': {'alpha':bnlearn_local_fdr, 'test':bnlearn_test, 'target_features':target_features, 
				'radius':max_radius, 'verbose':verbose}
		}

		method_args = bnlearn_local_methods[method_name]
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
	if mtype == 'huge':
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
	# SILGGM
	#----------------------------------------------------------------
	if mtype == 'silggm':
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
	if plot_result:
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
		filename = f"hcp_{method_name}_npn" if apply_npn else f"hcp_{method_name}"
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



