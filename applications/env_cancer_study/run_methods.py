# Apply PFS to county-level cancer data
"""
- Targets: Age-adjusted cancer incidence and mortality rates (all cancers)
- Covariates: environmental exposures, socioeconomic factors, and demographics
"""

# from graph_estimation.methods import run_method
from localgraph import max_cor_response, pfs, plot_graph, restrict_to_local_graph
import numpy as np
import pandas as pd
import pickle
from sklearn.preprocessing import StandardScaler

import sys, os
sys.path.insert(0, os.path.abspath('..'))
from methods import run_method

#----------------------------------------------------------------
# Setup
#----------------------------------------------------------------
save_result = True
plot_result = True

random_seed = 4161932
np.random.seed(random_seed)

methods_to_run = ['hpc_local']

apply_npn = False
max_cor_lambda = False
max_radius = 2
fdr_list = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1]
bnlearn_fdr = 0.05

#----------------------------------------------------------------
# Load cleaned data
#----------------------------------------------------------------
df = pd.read_csv(f'./data/cleaned_data/cleaned_data.csv')
X_raw = df.to_numpy()
feature_names = df.columns.tolist()
target_features = [feature_names.index(name) for name in ['Mortality', 'Incidence']]

#----------------------------------------------------------------
# Loop through methods
#----------------------------------------------------------------
for method_name in methods_to_run:

	print(f'Starting {method_name}')

	X = X_raw.copy()

	#----------------------------------------------------------------
	# bnlearn
	#----------------------------------------------------------------
	bnlearn_test = 'mi-g'
	bnlearn_methods = {
		'aracne':{'mi':'mi-g'},
		'fast_iamb':{'alpha':bnlearn_fdr, 'test':bnlearn_test},
		'hpc':{'alpha':bnlearn_fdr, 'test':bnlearn_test},
		'iamb':{'alpha':bnlearn_fdr, 'test':bnlearn_test},
		'iamb_fdr':{'alpha':bnlearn_fdr, 'test':bnlearn_test},
		'mmpc':{'alpha':bnlearn_fdr, 'test':bnlearn_test},
		'pc_stable':{'alpha':bnlearn_fdr, 'test':bnlearn_test},
		'si_hiton_pc':{'alpha':bnlearn_fdr, 'test':bnlearn_test}
	}
	if method_name in bnlearn_methods:
		method_args = bnlearn_methods[method_name]
		result = run_method(method_name, X, **method_args)
		A = result['adjacency_matrix']
		if max_radius is not None:
			result['adjacency_matrix'] = restrict_to_local_graph(A, target_features, max_radius)
		result['max_cor_lambda'] = max_cor_lambda
		runtime = result['runtime']
		print(f'Runtime: {runtime:.2f} seconds')

	#----------------------------------------------------------------
	# bnlearn local
	#----------------------------------------------------------------
	bnlearn_test = 'mi-g'
	bnlearn_local_methods = {
		'fast_iamb_local':{'alpha':bnlearn_fdr, 'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features, 'verbose':True},
		'hpc_local':{'alpha':bnlearn_fdr, 'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features, 'verbose':True},
		'iamb_local':{'alpha':bnlearn_fdr, 'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features, 'verbose':True},
		'mmpc_local':{'alpha':bnlearn_fdr, 'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features, 'verbose':True},
		'pc_stable_local':{'alpha':bnlearn_fdr, 'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features, 'verbose':True},
		'si_hiton_pc_local':{'alpha':bnlearn_fdr, 'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features, 'verbose':True}
	}
	if method_name in bnlearn_local_methods:
		method_args = bnlearn_local_methods[method_name]
		result = run_method(method_name, X, **method_args)
		A = result['adjacency_matrix']
		if max_radius is not None:
			result['adjacency_matrix'] = restrict_to_local_graph(A, target_features, max_radius)
		result['max_cor_lambda'] = max_cor_lambda
		runtime = result['runtime']
		print(f'Runtime: {runtime:.2f} seconds')

	#----------------------------------------------------------------
	# huge
	#----------------------------------------------------------------
	if method_name in ['glasso', 'mb']:
		X = StandardScaler().fit_transform(X)
		if max_cor_lambda:
			max_cors = max_cor_response(X, target_features)
			lambda_ = 0.99 * min(max_cors)
		else:
			lambda_ = None
		method_args = {'lambda_':lambda_, 'apply_npn':apply_npn}
		result = run_method(method_name, X, **method_args)
		A = result['adjacency_matrix']
		if max_radius is not None:
			result['adjacency_matrix'] = restrict_to_local_graph(A, target_features, max_radius)
		result['max_cor_lambda'] = max_cor_lambda
		runtime = result['runtime']
		print(f'Runtime: {runtime:.2f} seconds')

	#----------------------------------------------------------------
	# mixed graphical models (mgm)
	#----------------------------------------------------------------
	if method_name == 'mgm':
		X = StandardScaler().fit_transform(X)
		method_args = {'cat_threshold':4}
		result = run_method(method_name, X, **method_args)
		A = result['adjacency_matrix']
		if max_radius is not None:
			result['adjacency_matrix'] = restrict_to_local_graph(A, target_features, max_radius)
		result['max_cor_lambda'] = 'N/A'
		runtime = result['runtime']
		print(f'Runtime: {runtime:.2f} seconds')

	#----------------------------------------------------------------
	# SILGGM
	#----------------------------------------------------------------
	if method_name in ['bnwsl', 'dsgl', 'dsnwsl', 'gfcl', 'gfcsl']:
		X = StandardScaler().fit_transform(X)
		method_args = {'alpha':fdr_list, 'apply_npn':apply_npn}
		result = run_method(method_name, X, **method_args)
		adjacency_matrix = result['adjacency_matrix']
		target_fdrs = result['target_fdrs']
		if max_radius is not None:
			for target_fdr in target_fdrs:
				A = adjacency_matrix[target_fdr]
				result['adjacency_matrix'][target_fdr] = restrict_to_local_graph(A, target_features, max_radius)
		result['max_cor_lambda'] = 'N/A'
		runtime = result['runtime']
		print(f'Runtime: {runtime:.2f} seconds')

	#----------------------------------------------------------------
	# Plot result
	#----------------------------------------------------------------
	if plot_result:
		plot_args = {
			'target_features':target_features,
			'feature_names':feature_names,
			'radius':max_radius, 
			'font_size':10,
			'node_size':2000,
			'edge_widths':3,
			'show_weights':False
			}
		
		if method_name not in ['bnwsl', 'dsgl', 'dsnwsl', 'gfcl', 'gfcsl']:
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
		filename = f"eqi_{method_name}_npn.pkl" if apply_npn else f"eqi_{method_name}.pkl"
		result['method_args'] = method_args
		result['target_features'] = target_features
		result['feature_names'] = feature_names
		result['apply_npn'] = apply_npn
		with open(filename, "wb") as f:
			pickle.dump(result, f)
		print(f"Saved result to {filename}")

	print(f'Finished {method_name}')


















