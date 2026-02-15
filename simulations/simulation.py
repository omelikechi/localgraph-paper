# Compare different methods

import logging
import pickle
import time

from localgraph import pfs, tp_and_fp
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

import sys, os
sys.path.insert(0, os.path.abspath('..'))
from methods import run_method
from simulate_block import block_graph

#----------------------------------------------------------------
# Global settings
#----------------------------------------------------------------
################################
save_results = False
################################
do_dense = False
################################
do_nonlinear = False
################################
file_name = f'sim_results_nonlinear' if do_nonlinear else 'sim_results_linear'
file_name += '_dense' if do_dense else ''
################################
default_settings = False
################################

# random seed list
random_seed_list = np.arange(1,5)

#----------------------------------------------------------------
# Simulation parameters
#----------------------------------------------------------------

if default_settings:

	from default_settings import default_settings

	target_type = 'nonlinear' if do_nonlinear else 'linear'
	sparsity = 'dense' if do_dense else 'sparse'
	n = 100

	key = f'{target_type}_{sparsity}_n{n}'
	cfg = default_settings[key].copy()
	globals().update(cfg)

else:
	n = 100
	snr = 4 if do_nonlinear else None
	block_sizes = [1, 4, 195]
	p = sum(block_sizes)
	block_degree = [0, 0, 2]
	connector_degree = [4, 6]
	block_magnitude = np.ones(len(block_sizes))
	connector_magnitude = np.ones(len(block_sizes) - 1)
	lmin = 0.01
	lmax = 10
	radii = [1, 2]
	qpath_max = 0.2
	fdr_local = [0.2, 0.1, 0.1, 0.1]
	fdr = 0.1
	ipss_selector = 'gb' if do_nonlinear else 'adaptive_lasso' if n > p else 'l1'

# methods to run
silggm_methods = ['bnwsl', 'dsgl', 'dsnwsl', 'gfcl', 'gfcsl']
methods = ['inter_iamb_local']

# append method names to file name
file_name += f'_n{n}'
file_name = f'{file_name}_' + '_'.join(methods)

# simulation metadata, saved in the final output to refer back to
simulation_metadata = {'n':n, 'p':p, 'snr':snr, 'block_sizes':block_sizes, 'block_degree':block_degree, 'connector_degree':connector_degree,
'block_magnitude':block_magnitude, 'connector_magnitude':connector_magnitude, 'lmin':lmin, 'lmax':lmax, 
'random_seed_list':random_seed_list.tolist(), 'radii':radii, 'methods':methods, 'qpath_max':qpath_max, 'fdr_local':fdr_local,
'fdr':fdr, 'ipss_selector':ipss_selector, 'do_nonlinear':do_nonlinear}

# store results
all_results = []

# create nonlinear target if do_nonlinear is True
def nonlinear_target(X, A_true, target, neighbors, snr):
	signal = np.zeros(X.shape[0])
	for i in neighbors:
		signal += np.exp(-X[:,i]**2 / 2)
		A_true[target,i] = A_true[i,target] = 1
	sigma2 = np.var(signal) / snr
	X[:,target] = signal + np.random.normal(0, np.sqrt(sigma2), size=n) 
	# add edges between neighbors of target
	for i in neighbors:
		for j in neighbors:
			if i < j:
				A_true[i,j] = A_true[j,i] = 1
	return X, A_true

# convert dictionary of q-values to adjacency matrix
def dict_to_matrix(graph_dict, p):
	Q = graph_dict
	A = np.zeros((p,p))
	for (i,j), q in Q.items():
		A[i,j] = q
		A[j,i] = q
	return A

#----------------------------------------------------------------
# Method configurations
#----------------------------------------------------------------
lambda_ = None
huge_crit = 'ric' # options: ebic, stars, ric
bnlearn_test = 'mi-g' if do_nonlinear else 'cor'
criterion = 'forward'
verbose = False if save_results else True

method_configs = {
	# bnlearn (global)
	'aracne':{'mi':'mi-g'},
	'fast_iamb':{'alpha':fdr, 'test':bnlearn_test},
	'hpc':{'alpha':fdr, 'test':bnlearn_test},
	'iamb':{'alpha':fdr, 'test':bnlearn_test},
	'iamb_fdr':{'alpha':fdr, 'test':bnlearn_test},
	'inter_iamb_fdr':{'alpha':fdr, 'test':bnlearn_test},
	'mmpc':{'alpha':fdr, 'test':bnlearn_test},
	'pc_stable':{'alpha':fdr, 'test':bnlearn_test},
	'si_hiton_pc':{'alpha':fdr, 'test':bnlearn_test},

	# bnlearn (local)
	'fast_iamb_local':{'alpha':fdr, 'test':bnlearn_test, 'radius':max(radii), 'criterion':criterion, 'verbose':verbose},
	'hpc_local':{'alpha':fdr, 'test':bnlearn_test, 'radius':max(radii), 'criterion':criterion, 'verbose':verbose},
	'iamb_local':{'alpha':fdr, 'test':bnlearn_test, 'radius':max(radii), 'criterion':criterion, 'verbose':verbose},
	'iamb_fdr_local':{'alpha':fdr, 'test':bnlearn_test, 'radius':max(radii), 'criterion':criterion, 'verbose':verbose},
	'inter_iamb_local':{'alpha':fdr, 'test':bnlearn_test, 'radius':max(radii), 'criterion':criterion, 'verbose':verbose},
	'mmpc_local':{'alpha':fdr, 'test':bnlearn_test, 'radius':max(radii), 'criterion':criterion, 'verbose':verbose},
	'pc_stable_local':{'alpha':fdr, 'test':bnlearn_test, 'radius':max(radii), 'criterion':criterion, 'verbose':verbose},
	'si_hiton_pc_local':{'alpha':fdr, 'test':bnlearn_test, 'radius':max(radii), 'criterion':criterion, 'verbose':verbose},

	# huge
	'glasso':{'lambda_':lambda_, 'criterion':huge_crit},
	'mb':{'lambda_':None, 'criterion':huge_crit},

	# pfs
	'pfs':{'selector':ipss_selector, 'qpath_max':qpath_max, 'max_radius':max(radii), 
		'fdr_local':fdr_local, 'criterion':criterion, 'verbose':verbose},

	# silggm
	'bnwsl':{'alpha':fdr},
	'dsgl':{'alpha':fdr},
	'dsnwsl':{'alpha':fdr},
	'gfcl':{'alpha':fdr},
	'gfcsl':{'alpha':fdr},
}

simulation_metadata['method_configs'] = method_configs

#----------------------------------------------------------------
# Run simulation
#----------------------------------------------------------------
if save_results:
	logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s - %(message)s')

n_trials = len(random_seed_list)

print(f'Starting {file_name}')
print(f'----------------------------------------------------------------')

for trial, random_seed in enumerate(random_seed_list):

	if save_results:
		logging.info(f'trial {trial + 1}/{n_trials}')
	else:
		print(f'trial {trial + 1}/{n_trials}')

	# generate data
	data = block_graph(n, lmin, lmax, block_sizes, block_degree, block_magnitude, 
		connector_degree, connector_magnitude, random_seed=random_seed)

	# data output
	X = data['X']
	A_true = data['A']
	target_features = data['target_features']
	feature_names = data['feature_names']
	max_cor_response = data['max_cor_response']

	n_responses = len(target_features)
	n, p = X.shape

	# apply nonlinearity; note that in this study, the target feature is always 0
	if do_nonlinear:
		X, A_true = nonlinear_target(X, A_true, 0, np.arange(1, block_sizes[1]+1), snr)

	X = StandardScaler().fit_transform(X)		

	# add target_features to pfs and bnlearn_local args
	method_configs['pfs']['target_features'] = target_features
	for method_name in methods:
		if '_local' in method_name:
			method_configs[method_name]['target_features'] = target_features

	for i, method_name in enumerate(methods):
		if method_name == 'truth':
			continue
		else:
			config = method_configs[method_name]
			start = time.time()

			if method_name == 'pfs':
				Q = pfs(X, **config)
				A = dict_to_matrix(Q,p)
			else:
				result = run_method(method_name, X, **config)
				A = result['adjacency_matrix']

			A = np.maximum(A, A.T)
			method_time = time.time() - start

			print(f'  {method_name}: {method_time:.2f} seconds')

			# true positives in the full graph
			n_true_global = np.sum(A_true) // 2
			tp_global, fp_global = tp_and_fp(A, A_true, target_features, radius=None)
			tpr_global = tp_global / max(n_true_global, 1)
			fdp_global = fp_global / max(tp_global + fp_global, 1)

			# loop over radii
			for radius in radii:
				# count true edges eligible under the tp_and_fp definition
				tp_true, _ = tp_and_fp(A_true, A_true, target_features, radius=radius)
				n_true_local = tp_true

				tp_local, fp_local = tp_and_fp(A, A_true, target_features, radius=radius)
				tpr_local = tp_local / max(n_true_local, 1)
				fdp_local = fp_local / max(tp_local + fp_local, 1)

				# save per radius
				all_results.append({
					'seed': random_seed,
					'radius': radius,
					'method': method_name,
					'TPR_global': tpr_global,
					'FDP_global': fdp_global,
					'TPR_local': tpr_local,
					'FDP_local': fdp_local,
					'time_sec': method_time
				})

#----------------------------------------------------------------
# Save or print results summary
#----------------------------------------------------------------
# convert to dataframe and save
if save_results:
	results_package = {
		'metadata': simulation_metadata,
		'results': all_results
	}
	with open(f"{file_name}.pkl", "wb") as f:
		pickle.dump(results_package, f)
		logging.info("Simulation results saved.")
else:
	df_results = pd.DataFrame(all_results)
	radii = sorted(df_results['radius'].unique())
	methods = df_results['method'].unique()

	# Local TPR
	tpr_local_table = pd.DataFrame(index=radii, columns=methods)
	for r in radii:
		df_r = df_results[df_results['radius'] == r]
		for m in methods:
			vals = df_r[df_r['method'] == m]['TPR_local']
			mean = vals.mean()
			std = vals.std()
			tpr_local_table.loc[r, m] = f"{mean:.2f} ({std:.2f})"

	print("\nLocal TPR (mean (std))")
	print("----------------------------------------------------------------")
	print(tpr_local_table.to_string())

	# Local FDP
	fdp_local_table = pd.DataFrame(index=radii, columns=methods)
	for r in radii:
		df_r = df_results[df_results['radius'] == r]
		for m in methods:
			vals = df_r[df_r['method'] == m]['FDP_local']
			mean = vals.mean()
			std = vals.std()
			fdp_local_table.loc[r, m] = f"{mean:.2f} ({std:.2f})"

	print("\nLocal FDP (mean (std))")
	print("----------------------------------------------------------------")
	print(fdp_local_table.to_string())

print()


