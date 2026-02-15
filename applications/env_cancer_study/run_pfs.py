# Apply PFS to county-level cancer data
"""
- Targets: Age-adjusted cancer incidence and mortality rates (all cancers)
- Covariates: environmental exposures, socioeconomic factors, and demographics
- Note: PFS is a stochastic algorithm, so results may vary slightly across machines
	even when using the same random seed.
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
from localgraph import pfs, plot_graph

#--------------------------------
# Setup
#--------------------------------
save_result = False
show_result = True

random_seed = 4161932
np.random.seed(random_seed)

# Graph details
radius = 1
qpath_max = 1
fdr_local = [0.1, 0.008]

# custom neighborhood
"""
- Education: Many variables with the same q-value. Used efp scores to break ties. Custom variables ordered
	in terms of efp scores (Poverty: 0.04, Income: 0.04, HomeValue: 0.05, EntBusiness: 0.05, Colonoscopy: 0.06
	ViolentCrime: 0.06, White: 0.06, Mortality: 0.07)
"""
# custom_nbhd = {
# 	'Mortality':{'nbhd_fdr':0.009},
# 	'Education':{'nbhd_fdr':0.001, 'Poverty':0.1, 'Income':0.1, 'HomeValue':0.1, 'EntBusiness':0.1, 
# 		'Colonoscopy':0.1, 'ViolentCrime':0.1, 'White':0.1, 'Mortality':0.1},
# 	'Hispanic':{'nbhd_fdr':0.01},
# 	'pm2.5(A)':{'Education':0.001},
# 	'Poverty':{'nbhd_fdr':0.011},
# 	'Smoking':{'nbhd_fdr':0.01}
# }

custom_nbhd = {}

#--------------------------------
# Load cleaned data
#--------------------------------
data_dir = 'data/cleaned_data'
df = pd.read_csv(f'{data_dir}/cleaned_data.csv')
X = df.to_numpy()
feature_names = df.columns.tolist()
target_features = [feature_names.index(name) for name in ['Mortality', 'Incidence']]

#--------------------------------
# Run PFS
#--------------------------------
ipss_selector = 'rf'
Q = pfs(X, target_features, qpath_max=qpath_max, fdr_local=fdr_local, max_radius=radius, selector=ipss_selector, 
		custom_nbhd=custom_nbhd, feature_names=feature_names, verbose=True)

result = {'random_seed':random_seed, 'response_names':['Mortality',  'Incidence'], 'feature_names':feature_names, 
	'target_features':target_features, 'Q':Q}

# add metadata
result['X_shape'] = X.shape
result['custom_nbhd'] = custom_nbhd
result['fdr_local'] = fdr_local
result['radius'] = radius

#--------------------------------
# Save result
#--------------------------------
if save_result:
	output_path = f'results/env_cancer_pfs_results.pkl'
	with open(output_path, 'wb') as f:
		pickle.dump(result, f)

#--------------------------------
# Plot result
#--------------------------------
if show_result:
	fig, ax = plt.subplots(figsize=(17,6))

	plot_graph(
		result['Q'],
		ax=ax,
		target_features=target_features,
		feature_names=feature_names,
		radius=radius,
		node_size=3000,
		font_size=8,
		show_weights=True,
		graph_layout='kk',
		edge_digits=4
	)

	fig.tight_layout()
	plt.show()
