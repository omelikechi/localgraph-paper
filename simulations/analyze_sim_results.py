# Analyze simulation results and print LaTeX tables

import numpy as np
import pickle
import pandas as pd

"""
Methods that timed out:
  - SIHPC timed out (26 hours) in the n=500 linear dense simulation
  - SIHPC timed out (26 hours) in the n=500 nonlinear dense simulation
"""

################################
do_nonlinear = True
################################
do_dense = True
################################
n = 500
################################

simulation_type = 'nonlinear' if do_nonlinear else 'linear'
sparsity = 'dense' if do_dense else 'sparse'

if do_dense:
	file_name = f'./results/sim_results_{simulation_type}_dense_n{n}.pkl'
else:
	file_name = f'./results/sim_results_{simulation_type}_n{n}.pkl'

show_stdev = False
blue_text = True

#--------------------------------
# Load results
#--------------------------------
with open(file_name, "rb") as f:
	results_package = pickle.load(f)

metadata = results_package["metadata"]
for variable, value in metadata.items():
	if variable not in ['random_seed_list', 'methods']:
		print(f'{variable}: {value},')
print()
df = pd.DataFrame(results_package["results"])

df['method'] = df['method'].replace({'pfs': 'pfs', 'gipss': 'pfs'})

#--------------------------------
# Process
#--------------------------------
radii = sorted(df['radius'].unique())

# LaTeX method name mapping and order
latex_names = {
	'pfs':r'\textbf{PFS}',
	'glasso':'GLasso',
	'mb':'NLasso',
	'bnwsl':'BNWSL',
	'dsgl':'DSGL',
	'dsnwsl':'DSNWSL',
	'gfcl':'GFCL',
	'gfcsl':'GFCSL',
	'aracne':'ARACNE',
	'hpc':'HPC',
	'mmpc':'MMPC',
	'si_hiton_pc':'SIHPC',

	'fast_iamb_local':'FastIAMB(L)',
	'hpc_local':'HPC(L)',
	'iamb_local':'IAMB(L)',
	'mmpc_local':'MMPC(L)',
	'si_hiton_pc_local':'SIHPC(L)',
}
ordered_methods = ['pfs', 'glasso', 'mb', 'bnwsl', 'dsgl', 'dsnwsl', 'gfcl', 'gfcsl', 
	'aracne', 'hpc', 'mmpc', 'si_hiton_pc', 'hpc_local', 'mmpc_local', 'si_hiton_pc_local']

if not do_nonlinear and not do_dense and n == 100:
	ordered_methods += ['fast_iamb_local', 'iamb_local']

# process nan values due to method timing out after 24 hours
def fmt_mean(vals, show_stdev=False):
	mean = vals.mean()
	if np.isnan(mean):
		return r'--'
	if show_stdev:
		std = vals.std()
		return f'{mean:.2f} ({std:.2f})'
	return f'{mean:.2f}'

#--------------------------------
# Helper to format LaTeX table
#--------------------------------
def format_latex_table_tpr_fdr(df, caption, label, show_stdev=False):
	print("\\begin{table}[ht]")
	print("\\centering")
	if blue_text:
		print(r"{\color{blue}")
	print("\\begin{tabular}{l" + "c"*len(radii) + "c" + "c"*len(radii) + "}")
	print("\\toprule")
	print(" & \\multicolumn{" + str(len(radii)) + "}{c}{TPR} & & \\multicolumn{" + str(len(radii)) + "}{c}{FDR} \\\\")
	print("\\cmidrule(lr){2-" + str(1+len(radii)) + "}\\cmidrule(lr){" + str(3+len(radii)) + "-" + str(2+2*len(radii)) + "}")
	header = ["Method"] + [f"$r={r}$" for r in radii] + [""] + [f"$r={r}$" for r in radii]
	print(" & ".join(header) + " \\\\")
	print("\\midrule")
	for m in ordered_methods:
		row = [latex_names[m]]
		for r in radii:
			vals = df[(df['method'] == m) & (df['radius'] == r)]['TPR_local']
			row.append(fmt_mean(vals, show_stdev))
		row.append("")
		for r in radii:
			vals = df[(df['method'] == m) & (df['radius'] == r)]['FDP_local']
			row.append(fmt_mean(vals, show_stdev))
		print(" & ".join(row) + " \\\\")
	print("\\bottomrule")
	print("\\end{tabular}")
	if blue_text:
		print(r"}")
		print(f"\\caption{{\\blue{{\\textit{{{caption}}}")
	else:
		print(f"\\caption{{\\textit{{{caption}")
	print(f"\n%Simulation details: {metadata_str}.'")
	print(f"\\label{{{label}}}")
	print("\\end{table}")
	print()

#--------------------------------
# Print LaTeX table
#--------------------------------
metadata_str = ', '.join(
	f'{k}: {v}'
	for k, v in metadata.items()
	if k not in ['methods', 'random_seed_list']
)

format_latex_table_tpr_fdr(
	df=df,
	caption=(
		f'({simulation_type.capitalize()}, {sparsity}, $n = {n}$)}}. '
		f'True positive rate (TPR) and false discovery rate (FDR) '
		f'for local graph recovery across methods (rows) and neighborhood radii (columns), '
		f'averaged over 100 trials. Methods with a trailing (L) are local versions of '
		f'their global couterparts; for example, HPC(L) is the local version of HPC.}}'
	),
	label=f"tab:{simulation_type}_{sparsity}_n{n}",
	show_stdev=show_stdev
)




