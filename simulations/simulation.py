# Compare different methods

import logging
import pickle
import time

from localgraph import pfs, subgraph_within_radius, tp_and_fp
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

from methods import hugeR, silggmR
from simulate_block import block_graph

#--------------------------------
# Global settings
#--------------------------------
################################
save_results = False
################################
do_nonlinear = False
################################
file_name = f'simulation_results_nonlinear' if do_nonlinear else 'simulation_results_linear'
################################

# random seed list
random_seed_list = np.arange(1,4)

# simulation parameters
n = 100
lmin = 0.01
lmax = 10
block_sizes = [1, 4, 195]
block_degree = [0, 0, 2]
connector_degree = [4, 2]
block_magnitude = np.ones(len(block_sizes))
connector_magnitude = np.ones(len(block_sizes) - 1)
snr = 4

# number of variables is the sum of the block sizes
p = 0
for i in block_sizes:
	p += i

# radii at which to evaluate local graph estimation performance
radii = [1,2]

# methods to run
methods = ['truth', 'glasso', 'mb', 'dsnwsl', 'dsgl', 'bnwsl', 'gfcsl', 'gfcl', 'pfs']

# ipss selector for pfs ('l1': use lasso as base estimator; 'gb': use gradient boosting)
ipss_selector = 'gb' if do_nonlinear else 'l1'

# false discovery parameters for pfs
qpath_max = 0.2
fdr_local = [0.2, 0.05, 0.05, 0.05] if do_nonlinear else [0.2, 0.1, 0.1, 0.1]

# target FDR for other methods
fdr = 0.1

# simulation metadata, saved in the final output to refer back to
simulation_metadata = {'n':n, 'p':p, 'snr':snr, 'block_sizes':block_sizes, 'block_degree':block_degree, 'connector_degree':connector_degree,
'block_magnitude':block_magnitude, 'connector_magnitude':connector_magnitude, 'lmin':lmin, 'lmax':lmax, 
'random_seed_list':random_seed_list.tolist(), 'radii':radii, 'methods':methods, 'qpath_max':qpath_max, 'fdr_local':fdr_local,
'fdr':fdr, 'ipss_selector':ipss_selector, 'do_nonlinear':do_nonlinear}

# store results
all_results = []

# create a nonlinear target if do_nonlinear is True
def nonlinear_target(X, A_true, target, neighbors, snr):
	signal = np.zeros(X.shape[0])
	for i in neighbors:
		signal += np.exp(-X[:,i]**2 / 2)
		A_true[target,i] = A_true[i,target] = 1
	sigma2 = np.var(signal) / snr
	X[:,target] = signal + np.random.normal(0, np.sqrt(sigma2), size=n) 
	return X, A_true

# convert dictionary of q-values to adjacency matrix
def dict_to_matrix(graph_dict, p):
	Q = graph_dict
	A = np.zeros((p,p))
	for (i,j), q in Q.items():
		A[i,j] = q
		A[j,i] = q
	return A

#--------------------------------
# Method configurations
#--------------------------------
lambda_ = None
flare_crit = 'stars' # options: cv, stars
huge_crit = 'ric' # options: ebic, stars, ric
method_configs = {
	'bnwsl': {'method':'B_NW_SL', 'alpha':fdr},
	'dsgl':{'method':'D-S_GL', 'alpha':fdr},
	'dsnwsl':{'method':'D-S_NW_SL', 'alpha':fdr},
	'gfcl':{'method':'GFC_L', 'alpha':fdr},
	'gfcsl':{'method':'GFC_SL', 'alpha':fdr},
	'glasso':{'method':'glasso', 'sym':'or', 'lambda_':lambda_, 'criterion':huge_crit},
	'mb':{'method':'mb', 'sym':'or', 'lambda_':None, 'criterion':huge_crit},
	'pfs':{'selector':ipss_selector, 'qpath_max':qpath_max, 'max_radius':max(radii), 
		'fdr_local':fdr_local, 'criterion':'forward'}
}

# function that runs methods
def run_method(X, method_name, **kwargs):
	method_functions = {
		'bnwsl': lambda X, kwargs: silggmR(X, **kwargs),
		'dsgl': lambda X, kwargs: silggmR(X, **kwargs),
		'dsnwsl': lambda X, kwargs: silggmR(X, **kwargs),
		'gfcl': lambda X, kwargs: silggmR(X, **kwargs),
		'gfcsl': lambda X, kwargs: silggmR(X, **kwargs),
		'glasso': lambda X, kwargs: hugeR(X, **kwargs),
		'mb': lambda X, kwargs: hugeR(X, **kwargs),
		'pfs': lambda X, kwargs: pfs(X, **kwargs)
	}
	if method_name == 'pfs':
		Q = method_functions[method_name](X, kwargs)
		A = dict_to_matrix(Q,p)
		return A
	else:
		return method_functions[method_name](X, kwargs)

#--------------------------------
# Run simulation
#--------------------------------
if save_results:
	logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s - %(message)s')

n_trials = len(random_seed_list)
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
		X, A_true = nonlinear_target(X, A_true, 0, np.arange(1,block_sizes[1]+1), snr)

	X = StandardScaler().fit_transform(X)

	# add target_features to pfs args
	method_configs['pfs']['target_features'] = target_features

	for i, method_name in enumerate(methods):
		if method_name == 'truth':
			continue
		else:
			config = method_configs[method_name]
			start = time.time()
			A = run_method(X, method_name, **config)
			A = np.maximum(A, A.T)
			method_time = time.time() - start

			# true positives in the full graph
			n_true_global = np.sum(A_true) // 2
			tp_global, fp_global = tp_and_fp(A, A_true, target_features, radius=None)
			tpr_global = tp_global / max(n_true_global, 1)
			fdp_global = fp_global / max(tp_global + fp_global, 1)

			# loop over radii
			for radius in radii:
				A_true_local = subgraph_within_radius(A_true, target_features, radius)
				n_true_local = np.sum(A_true_local) // 2

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

#--------------------------------
# Save or print results summary
#--------------------------------
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
	print("------------------------------------------------------------------------------------------------")
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
	print("------------------------------------------------------------------------------------------------")
	print(fdp_local_table.to_string())
