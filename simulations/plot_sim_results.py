import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

################################
save_fig = True
################################

#----------------------------------------------------------------
# Settings
#----------------------------------------------------------------
old_results = False
do_nonlinear = False
simulation_type = 'nonlinear' if do_nonlinear else 'linear'

n_list = [50, 100, 150, 200, 250, 300, 350, 400]
methods_to_plot = ['pfs'] #, 'glasso', 'mb']
results_dir = './results/old' if old_results else './results/'

#----------------------------------------------------------------
# Load and combine results
#----------------------------------------------------------------
dfs = []

for n in n_list:
	if n == 201:
		file_name = f'{results_dir}/sim_results_{simulation_type}_n{n}_pfs_adaptivelasso.pkl'
	else:
		file_name = f'{results_dir}/sim_results_{simulation_type}_n{n}.pkl'
	with open(file_name, 'rb') as f:
		output = pickle.load(f)

	df = pd.DataFrame(output['results'])
	df['n'] = n
	df['method'] = df['method'].replace({'gipss': 'pfs'})

	print(f'n = {n}, ipss_selector = {output['metadata']['ipss_selector']}')

	dfs.append(df)

df_all = pd.concat(dfs, ignore_index=True)

# keep only selected methods
df_all = df_all[df_all['method'].isin(methods_to_plot)]

radii = sorted(df_all['radius'].unique())

# average over trials
summary = (
	df_all
	.groupby(['method', 'radius', 'n'])
	.agg(
		TPR=('TPR_local', 'mean'),
		FDR=('FDP_local', 'mean')
	)
	.reset_index()
)

#----------------------------------------------------------------
# Plot
#----------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(16,6), sharex=True, sharey=True)

# fdr vs n
for m in methods_to_plot:
	for r in radii:
		df_mr = summary[(summary['method'] == m) & (summary['radius'] == r)]
		axes[0].plot(
			df_mr['n'],
			df_mr['FDR'],
			marker='o',
			label=f'{m}, r={r}',
			lw=3,
			markersize=8
		)

axes[0].set_xlabel(f'$n$', fontsize=28)
axes[0].set_ylabel('FDR', fontsize=28)

# tpr vs n
for m in methods_to_plot:
	for r in radii:
		df_mr = summary[(summary['method'] == m) & (summary['radius'] == r)]
		axes[1].plot(
			df_mr['n'],
			df_mr['TPR'],
			marker='o',
			label=f'r={r}',
			lw=3,
			markersize=8
		)

axes[1].set_xlabel(f'$n$', fontsize=28)
axes[1].set_ylabel('TPR', fontsize=28)
axes[1].legend(loc='best', fontsize=22)

for ax in axes:
	ax.tick_params(axis='x', labelsize=16)
	ax.tick_params(axis='y', labelsize=12)

plt.tight_layout()
if save_fig:
	plt.savefig('varying_n.png', dpi=300)
plt.show()
