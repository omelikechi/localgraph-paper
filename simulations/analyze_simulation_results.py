# Analyze simulation results and print LaTeX tables

import pickle
import pandas as pd

# do_nonlinear = True for nonlinear results, False for linear results
do_nonlinear = False
simulation_type = 'nonlinear' if do_nonlinear else 'linear'
file_name = f'./results/simulation_results_{simulation_type}.pkl'

#--------------------------------
# Load results
#--------------------------------
with open(file_name, "rb") as f:
	results_package = pickle.load(f)

metadata = results_package["metadata"]
df = pd.DataFrame(results_package["results"])

#--------------------------------
# Process
#--------------------------------
radii = sorted(df['radius'].unique())

# LaTeX method name mapping and order
latex_names = {
	'gipss': r'\textbf{PFS}',
	'glasso': 'GLasso',
	'mb': 'NLasso',
	'bnwsl': 'BNWSL',
	'dsgl': 'DSGL',
	'dsnwsl': 'DSNWSL',
	'gfcl': 'GFCL',
	'gfcsl': 'GFCSL'
}
ordered_methods = ['gipss', 'glasso', 'mb', 'bnwsl', 'dsgl', 'dsnwsl', 'gfcl', 'gfcsl']

#--------------------------------
# Helper to format LaTeX table
#--------------------------------
def format_latex_table(df_metric, metric_name, caption, label):
	print("\\begin{table}[ht]")
	print("\\centering")
	print("\\begin{tabular}{ccccccccc}")
	print("\\toprule")
	header = ["Radius"] + [latex_names[m] for m in ordered_methods]
	print(" & ".join(header) + " \\\\")
	print("\\midrule")
	for r in radii:
		row = [str(r)]
		for m in ordered_methods:
			vals = df[(df['radius'] == r) & (df['method'] == m)][metric_name]
			mean = vals.mean()
			row.append(f"{mean:.2f}")
		print(" & ".join(row) + " \\\\")
	print("\\bottomrule")
	print("\\end{tabular}")
	print(f"\\caption{{\\textit{{{caption}}}}}")
	print(f"\\label{{{label}}}")
	print("\\end{table}")
	print()

#--------------------------------
# Print LaTeX tables
#--------------------------------
format_latex_table(
	df_metric=df,
	metric_name="TPR_local",
	caption=f"TPR by radius, {simulation_type} simulation study (average over 100 trials).",
	label=f"tab:sim_{simulation_type}_tpr"
)

format_latex_table(
	df_metric=df,
	metric_name="FDP_local",
	caption=f"FDR by radius, {simulation_type} simulation study (average over 100 trials).",
	label=f"tab:sim_{simulation_type}_fdr"
)


