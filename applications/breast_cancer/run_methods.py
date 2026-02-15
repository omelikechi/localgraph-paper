# Apply graph estimation methods to TCGA breast cancer data
"""
- Targets: histological type, pathologic stage, and status (alive or deceased)
- Covariates: genes (RNAseq), proteins (RPPA), and microRNAs
"""

from localgraph import max_cor_response, plot_graph, restrict_to_local_graph
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
save_result = False
plot_result = True

random_seed = 11201959
np.random.seed(random_seed)

drop_feature_types = ['rnaseq', 'mirna']

methods_to_run = ['fast_iamb_local']

apply_npn = False
max_cor_lambda = False
max_radius = 2
target_fdrs = [0.5] #[0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]
bnlearn_fdr = 0.05

huge_methods = ['glasso', 'mb']
silggm_methods = ['bnwsl', 'dsgl', 'dsnwsl', 'gfcl', 'gfcsl']

#----------------------------------------------------------------
# Load cleaned data
#----------------------------------------------------------------
df = pd.read_csv(f'./data/cleaned_data/cleaned_data.csv')

if drop_feature_types:
	mask = np.ones(len(df.columns), dtype=bool)
	for ftype in drop_feature_types:
		if ftype == 'rnaseq':
			mask &= ~df.columns.str.endswith('*')
		elif ftype == 'rppa':
			mask &= ~df.columns.str.endswith('#')
		elif ftype == 'mirna':
			mask &= ~(
				df.columns.str.contains('mir-', case=False, regex=False)
				| df.columns.str.contains('let-', case=False, regex=False)
			)
	df = df.loc[:, mask]

# Keep original column names
feature_names = df.columns.tolist()

X_raw = df.to_numpy()
target_names = ['histological_type', 'pathologic_stage', 'status']
target_features = [df.columns.get_loc(name) for name in target_names]

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
		'fast_iamb_local':{'alpha':bnlearn_fdr, 'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features},
		'hpc_local':{'alpha':bnlearn_fdr, 'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features},
		'iamb_local':{'alpha':bnlearn_fdr, 'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features},
		'mmpc_local':{'alpha':bnlearn_fdr, 'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features},
		'pc_stable_local':{'alpha':bnlearn_fdr, 'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features},
		'si_hiton_pc_local':{'alpha':bnlearn_fdr, 'test':bnlearn_test, 'radius':max_radius, 'target_features':target_features}
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
	if method_name in huge_methods:
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
	# Mixed graphical models (mgm)
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
	if method_name in silggm_methods:
		X = StandardScaler().fit_transform(X)
		method_args = {'alpha': target_fdrs, 'apply_npn': apply_npn}
		result = run_method(method_name, X, **method_args)
		adjacency_matrix = result['adjacency_matrix']
		if max_radius is not None:
			if len(target_fdrs) == 1:
				alpha = target_fdrs[0]
				A = adjacency_matrix
				result['adjacency_matrix'] = {alpha: restrict_to_local_graph(A, target_features, max_radius)}
			else:
				for alpha in target_fdrs:
					A = adjacency_matrix[alpha]
					result['adjacency_matrix'][alpha] = restrict_to_local_graph(A, target_features, max_radius)
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
		filename = f"bc_{method_name}_npn" if apply_npn else f"bc_{method_name}"
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


