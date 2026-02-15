# Analyze simulation results and append new methods (safe version)

import pickle
import pandas as pd

#--------------------------------
# Settings
#--------------------------------
do_nonlinear = True
do_dense = True
n = 500

simulation_type = 'nonlinear' if do_nonlinear else 'linear'

new_methods = ['hpc_local_mmpc_local_si_hiton_pc_local']

if do_dense:
	old_file = f'./results/sim_results_{simulation_type}_dense_n{n}.pkl'
	out_file = f'sim_results_{simulation_type}_dense_n{n}.pkl'
else:
	old_file = f'./results/sim_results_{simulation_type}_n{n}.pkl'
	out_file = f'sim_results_{simulation_type}_n{n}.pkl'

#--------------------------------
# Load existing results
#--------------------------------
with open(old_file, 'rb') as f:
	old_package = pickle.load(f)

metadata = old_package['metadata']
df_old = pd.DataFrame(old_package['results'])

dfs_new = []

#--------------------------------
# Load each new method package
#--------------------------------
for method in new_methods:
	if do_dense:
		file = f'./results/sim_results_{simulation_type}_dense_n{n}_{method}.pkl'
	else:
		file = f'./results/sim_results_{simulation_type}_n{n}_{method}.pkl'

	with open(file, 'rb') as f:
		pkg = pickle.load(f)

	df_new = pd.DataFrame(pkg['results'])

	assert set(df_old.columns) == set(df_new.columns)

	dfs_new.append(df_new)

#--------------------------------
# Append all at once
#--------------------------------
df_combined = pd.concat([df_old] + dfs_new, ignore_index=True)

#--------------------------------
# Update metadata
#--------------------------------
metadata = metadata.copy()
metadata['methods'] = sorted(df_combined['method'].unique())

#--------------------------------
# Save
#--------------------------------
results_package_combined = {
	'metadata': metadata,
	'results': df_combined.to_dict(orient='records')
}

with open(out_file, 'wb') as f:
	pickle.dump(results_package_combined, f)


