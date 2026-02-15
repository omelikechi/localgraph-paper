# Combine results for different dimensions into one file

import glob
import pandas as pd
import ast
import pickle

max_time = 2 * 60 * 60
p_list = [125, 250, 500, 1000, 2000, 4000, 8000]

files = glob.glob('./runtimes_by_dimension/runtime_test_*.csv')
if not files:
	raise RuntimeError('No runtime_test_*.csv files found.')

df_all = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)

def unwrap(x):
	if isinstance(x, str) and x.startswith('{'):
		d = ast.literal_eval(x)
		return list(d.values())[0]
	return x

for col in ['n', 'p_over_n', 'method', 'time_sec']:
	if col in df_all.columns:
		df_all[col] = df_all[col].apply(unwrap)

df_all['p'] = pd.to_numeric(df_all['p'], errors='coerce')
df_all['n'] = pd.to_numeric(df_all['n'], errors='coerce')
df_all['time_sec'] = pd.to_numeric(df_all['time_sec'], errors='coerce')

df_all = df_all.dropna(subset=['p', 'method', 'time_sec'])
df_all = df_all.sort_values(['method', 'p']).drop_duplicates(['method', 'p'], keep='last')

for method, g in df_all.groupby('method'):

	g = g.set_index('p')

	runtime_dict = {}
	for p in p_list:
		if p in g.index:
			runtime_dict[int(p)] = float(g.loc[p]['time_sec'])
		else:
			runtime_dict[int(p)] = float(max_time)

	n_val = int(g['n'].dropna().iloc[0]) if g['n'].notna().any() else None

	out = {
		'method': method,
		'n': n_val,
		'runtimes': runtime_dict
	}

	with open(f'runtime_{method}.pkl', 'wb') as f:
		pickle.dump(out, f)

	# print(out)


















# import glob
# import os
# import pandas as pd
# import ast

# max_time = 2 * 60 * 60
# p_list = [125, 250, 500, 1000, 2000, 4000, 8000]

# files = glob.glob('./runtimes_by_dimension/runtime_test_*.csv')
# if not files:
# 	raise RuntimeError('No runtime_test_*.csv files found.')

# df_all = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)

# print(df_all)

# def unwrap(x):
# 	if isinstance(x, str) and x.startswith('{'):
# 		d = ast.literal_eval(x)
# 		return list(d.values())[0]
# 	return x

# for col in ['n', 'p_over_n', 'method', 'time_sec']:
# 	if col in df_all.columns:
# 		df_all[col] = df_all[col].apply(unwrap)

# # **critical fixes**
# df_all['p'] = pd.to_numeric(df_all['p'], errors='coerce')
# df_all['n'] = pd.to_numeric(df_all['n'], errors='coerce')
# df_all['p_over_n'] = pd.to_numeric(df_all['p_over_n'], errors='coerce')
# df_all['time_sec'] = pd.to_numeric(df_all['time_sec'], errors='coerce')

# df_all = df_all.dropna(subset=['p', 'method', 'time_sec'])
# df_all = df_all.sort_values(['method', 'p']).drop_duplicates(['method', 'p'], keep='last')

# for method, g in df_all.groupby('method'):

# 	g = g.set_index('p')

# 	rows = []
# 	for p in p_list:
# 		if p in g.index:
# 			row = g.loc[p].to_dict()
# 			row['p'] = p
# 		else:
# 			row = {'method':method, 'p':p, 'time_sec':max_time}

# 		rows.append(row)

# 	df_method = pd.DataFrame(rows).sort_values('p')
# 	print(df_method)
# 	df_method.to_csv(f'runtime_{method}.csv', index=False)
















# import glob
# import os
# import pandas as pd

# max_time = 2 * 60 * 60
# p_list = [125, 250, 500, 1000, 2000, 4000, 8000]

# #----------------------------------------------------------------
# # Load all existing CSVs
# #----------------------------------------------------------------
# files = glob.glob('./runtimes_by_dimension/runtime_test_*.csv')
# if not files:
# 	raise RuntimeError('No runtime_test_*.csv files found.')

# df_all = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)

# #----------------------------------------------------------------
# # Aggregate per method and fill timeouts
# #----------------------------------------------------------------
# for method, g in df_all.groupby('method'):

# 	g = g.set_index('p')

# 	rows = []
# 	for p in p_list:
# 		if p in g.index:
# 			row = g.loc[p].to_dict()
# 			row['p'] = p
# 		else:
# 			row = {'method':method, 'p':p, 'time_sec':max_time}

# 		rows.append(row)

# 	df_method = pd.DataFrame(rows).sort_values('p')

# 	print(df_method)

# 	df_method.to_csv(f'runtime_{method}.csv', index=False)



