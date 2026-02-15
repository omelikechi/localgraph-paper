# Analyze method runtimes

import os
import pickle

import glob
import matplotlib.pyplot as plt
import pandas as pd

#----------------------------------------------------------------
# Setup
#----------------------------------------------------------------

show_plot = False

methods_to_plot = {
	'pfs_l1_r3':'\\textbf{PFS}(L1)',
	'aracne':'ARACNE',
	'pfs_gb_r3':'\\textbf{PFS}(GB)',
	'mb':'NLasso',
	'gfcsl':'GFCSL',
	'dsnwsl':'DSNWSL',
	'glasso':'GLasso',
	'mmpc':'MMPC',
	'si_hiton_pc':'SIHPC',
	'hpc':'HPC',
	'pc_stable':'PC',
	'fast_iamb':'FastIAMB',
	'hpc_local':'HPC(L)',
	'mmpc_local':'MMPC(L)',
	'si_hiton_pc_local':'SIHPC(L)',
	'fast_iamb_local':'FastIAMB(L)',
	'iamb_local':'IAMB(L)'
}

use_log_scale = False

files = glob.glob('./runtime_results/runtime_*.pkl')

#----------------------------------------------------------------
# Load all data
#----------------------------------------------------------------
data_by_method = {}

for f in files:
	with open(f, 'rb') as fp:
		obj = pickle.load(fp)

	method = obj['method']
	n = obj['n']
	runtimes = obj['runtimes']

	df = pd.DataFrame({
		'p': list(runtimes.keys()),
		'time_sec': list(runtimes.values())
	})

	data_by_method[method] = df

#----------------------------------------------------------------
# Plot runtime vs p
#----------------------------------------------------------------
plt.figure(figsize=(14,6))
all_rows = []

for method in methods_to_plot:
	if method not in data_by_method:
		continue

	df_sorted = data_by_method[method].sort_values('p')

	for _, row in df_sorted.iterrows():
		all_rows.append({'method': method, 'p': row['p'], 'time_sec': row['time_sec']})

	plt.plot(df_sorted['p'], df_sorted['time_sec'], marker='o', label=methods_to_plot[method])
	
if show_plot:
	plt.xlabel('p')
	plt.ylabel('Runtime (seconds)')

	if use_log_scale:
		plt.xscale('log')
		plt.yscale('log')

	plt.legend()
	plt.tight_layout()
	plt.show()

#----------------------------------------------------------------
# Create LaTeX table
#----------------------------------------------------------------
df_table = pd.DataFrame(all_rows)

def sec_to_hms(s):
	if s >= 2 * 60 * 60:
		return r'$>$2hrs'
	h = int(s // 3600)
	m = int((s % 3600) // 60)
	sec = int(s % 60)
	return f'{h}:{m:02d}:{sec:02d}'

table = df_table.pivot(index='method', columns='p', values='time_sec')

# sort columns
p_display = sorted(table.columns)
p_sort = sorted(table.columns, reverse=True)

table = table.reindex(columns=p_display)
table = table.sort_values(by=p_sort)

table_fmt = table.map(sec_to_hms)

# print table
p_cols = table_fmt.columns.tolist()

print(r'\begin{table*}[ht]')
print(r'\centering')
print(r'\begin{tabular}{l' + 'c' * len(p_cols) + r'}')
print(r'\toprule')

header = r'Method & ' + ' & '.join([f'{int(p)}' for p in p_cols]) + r' \\'
print(header)
print(r'\midrule')

for method in table_fmt.index:
	label = methods_to_plot[method]
	row = ' & '.join(table_fmt.loc[method].tolist())
	print(f'{label} & {row} \\\\')

print(r'\bottomrule')
print(r'\end{tabular}')
print(r'\end{table*}')

print()


