# Runtimes of different methods

import time

from localgraph import pfs
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

from default_settings import default_settings

import sys, os
sys.path.insert(0, os.path.abspath('../..'))
from methods import run_method
from methods.metadata import method_type
from simulations.simulate_block import block_graph

################################
save_results = False
################################

#----------------------------------------------------------------
# Global settings
#----------------------------------------------------------------
random_seed = 302
p_list = [125, 250, 500, 1000, 2000, 4000, 8000]
method = 'mmpc_local'

m_type = method_type(method)

# default or custom settings
use_default_settings = True
if use_default_settings:
	n = default_settings['n']
	method_args = default_settings[m_type + '_args']
else:
	n = 500
	method_args = {}

results = []

#----------------------------------------------------------------
# Run test
#----------------------------------------------------------------
for p in p_list:
	print(f'\n{method} (p = {p})')
	print(f'--------------------------------')

	# simple block construction: 1 target + rest noise
	block_sizes = [1, p - 1]
	block_degree = [0, 3]
	connector_degree = [3]
	block_magnitude = np.ones(len(block_sizes))
	connector_magnitude = np.ones(len(block_sizes) - 1)

	data = block_graph(
		n=n,
		lmin=0.01,
		lmax=10,
		block_sizes=block_sizes,
		block_degree=block_degree,
		block_magnitude=block_magnitude,
		connector_degree=connector_degree,
		connector_magnitude=connector_magnitude,
		random_seed=random_seed
	)

	X = StandardScaler().fit_transform(data['X'])
	target_features = data['target_features']

	if 'pfs' in method or m_type == 'bnlearn_local':
		method_args['target_features'] = target_features

	if 'pfs' in method:
		start = time.time()
		result = pfs(X, **method_args)
		runtime = time.time() - start
	else:
		result = run_method(method, X, **method_args)
		runtime = result['runtime']

	print(f'{method}: {runtime:.2f} seconds')

	results.append({'p':p, 'n':n, 'p_over_n':p / n, 'method':method, 'time_sec':runtime})

	#----------------------------------------------------------------
	# Save results
	#----------------------------------------------------------------
	if save_results:
		df = pd.DataFrame(results)
		df.to_csv(f'runtime_test_{method}_p{p}.csv', index=False)
		print(f'\nSaved runtime_test_{method}_p{p}.csv')


		
