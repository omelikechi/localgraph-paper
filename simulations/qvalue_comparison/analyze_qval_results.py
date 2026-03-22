# Analyze simulation results and print LaTeX tables

import numpy as np
import pickle
import pandas as pd

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
	file_name = f'../results/sim_results_{simulation_type}_dense_n{n}.pkl'
else:
	file_name = f'../results/sim_results_{simulation_type}_n{n}.pkl'

show_stdev = False
blue_text = True

#--------------------------------
# Load results
#--------------------------------
with open(file_name, "rb") as f:
	results_package = pickle.load(f)

metadata = results_package["metadata"]
# for variable, value in metadata.items():
# 	if variable not in ['random_seed_list', 'methods']:
# 		print(f'{variable}: {value},')
# print()
df = pd.DataFrame(results_package["results"])

df['method'] = df['method'].replace({'pfs': 'pfs', 'gipss': 'pfs'})

# print(metadata['methods'])

#--------------------------------
# Process
#--------------------------------
radii = sorted(df['radius'].unique())

# LaTeX method name mapping and order
latex_names = {
	'pfs':r'\textbf{IPSS}',
	'pfs_ipss_gb':r'\textbf{PFS(GB)}',
	'pfs_ipss_l1':r'\textbf{PFS(L1)}',
	'pfs_ipss_rf':r'\textbf{PFS(RF)}',
	'pfs_koglm':'KOGLM',
	'pfs_kol1':'KOL1',
	'pfs_korf':'KORF',
	'pfs_bh':'BH',
	'pfs_by':'BY',
}
ordered_methods = ['pfs', 'pfs_koglm', 'pfs_kol1', 'pfs_korf', 'pfs_bh', 'pfs_by']

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
		f'PFS with different $q$-value methods ({simulation_type.capitalize()}, {sparsity}, $n = {n}$)}}. '
		f'True positive rate (TPR) and false discovery rate (FDR) '
		f'for local graph recovery across methods (rows) and neighborhood radii (columns), '
		f'averaged over 100 trials. PFS with integrated path stability selection~\\cite{{ipss,ipss_nonparametric}} (IPSS) '
		f'is the method used throughout this work and the default option in the \\texttt{{localgraph}} software package. '
		f'Other methods are: model-X knockoffs~\\cite{{modelX}} with importance statistics based on (i) generalized linear ' 
		f'models (KOGLM), (ii) the lasso (KOL1), and (iii) random forests (KORF); the Benjamini--Hochberg method~\\cite{{bh}} ' 
		f'(BH); the Benjamini--Yekutieli method~\\cite{{by}} (BY).}}'
	),
	label=f"tab:qval_{simulation_type}_{sparsity}_n{n}",
	show_stdev=show_stdev
)




