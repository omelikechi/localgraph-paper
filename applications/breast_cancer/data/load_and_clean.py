# Functions for loading LinkedOmics data

import os
import numpy as np
import pandas as pd

pd.set_option('future.no_silent_downcasting', True)

def load_and_clean( 
	feature_types, 
	responses=None,
	remove_nan=True, 
	correlation_threshold=0.999, 
	expression_threshold=None,
	variance_threshold=None,
	print_feature_names=False, 
	verbose=False
):
	if verbose:
		print(f'Loading data')
		print(f'--------------------------------')

	base_path = os.path.dirname(os.path.dirname(__file__))
	file_paths = {
		'clinical': os.path.join(base_path, './data/raw_data/clinical.txt'),
		'mirna': os.path.join(base_path, './data/raw_data/mirna.txt'),
		'rnaseq': os.path.join(base_path, './data/raw_data/rnaseq.txt'),
		'rppa': os.path.join(base_path, './data/raw_data/rppa.txt')
	}

	shared_labels = None

	# Process responses
	if responses is not None:
		Y_df_list, shared_labels = process_responses(responses, file_paths, verbose=verbose)

	if print_feature_names:
		for category in feature_types:
			df = load_dataframe(file_paths[category])
			print(f'{category} features: {df.columns.tolist()}')
		return  # skip the rest if just printing features

	feature_dfs = []
	feature_names = []
	category_feature_names = {}

	for category in feature_types:
		df = load_dataframe(file_paths[category])

		if verbose:
			print(f' - Original {category} dimensions: {df.shape}')

		# Remove response columns if they exist
		if responses is not None:
			for response_type, response_name in responses:
				if category == response_type:
					df = df.drop(columns=[response_name], errors='ignore')

		# Update shared labels
		if shared_labels is not None:
			shared_labels = shared_labels.intersection(df.index)
		else:
			shared_labels = set(df.index)

		# Apply mean/variance filter on rnaseq
		# Apply mean/variance filter on rnaseq
		if category == 'rnaseq':

			mean_expression = df.mean(axis=0)
			variances = df.var(axis=0)

			if expression_threshold is None and variance_threshold is not None:
				var_cutoff = np.percentile(variances, variance_threshold)
				df = df.loc[:, variances > var_cutoff]
				if verbose:
					print(f' - {category} after variance filtering: {df.shape}')
					print(f'   - variance cutoff = {var_cutoff:.4f}')

			elif expression_threshold is not None and variance_threshold is None:
				mean_cutoff = np.percentile(mean_expression, expression_threshold)
				df = df.loc[:, mean_expression > mean_cutoff]
				if verbose:
					print(f' - {category} after mean filtering: {df.shape}')
					print(f'   - mean cutoff = {mean_cutoff:.2f}')

			elif expression_threshold is not None and variance_threshold is not None:
				mean_cutoff = np.percentile(mean_expression, expression_threshold)
				var_cutoff = np.percentile(variances, variance_threshold)
				joint_mask = (mean_expression > mean_cutoff) | (variances > var_cutoff)
				df = df.loc[:, joint_mask]
				if verbose:
					print(f' - {category} after joint mean/variance filtering: {df.shape}')
					print(f'   - mean cutoff = {mean_cutoff:.2f}, variance cutoff = {var_cutoff:.4f}')







		
		# if category == 'rnaseq':
		# 	if expression_threshold is None:

		# 	elif expression_threshold is not None or variance_threshold is not None:
		# 		expression_threshold = 0 if expression_threshold is None else expression_threshold
		# 		variance_threshold = 0 if variance_threshold is None else variance_threshold
		# 		mean_expression = df.mean(axis=0)
		# 		variances = df.var(axis=0)
		# 		mean_cutoff = np.percentile(mean_expression, expression_threshold)
		# 		var_cutoff = np.percentile(variances, variance_threshold)
		# 		joint_mask = (mean_expression > mean_cutoff) | (variances > var_cutoff)
		# 		df = df.loc[:, joint_mask]
		# 		if verbose:
		# 			print(f' - {category} after joint mean/variance filtering: {df.shape}')
		# 			print(f'   - mean cutoff = {mean_cutoff:.2f}, variance cutoff = {var_cutoff:.4f}')

		renamed_features = standardize_feature_names(category, df.columns.tolist(), feature_types)
		category_feature_names[category] = renamed_features
		feature_names.extend(renamed_features)
		df.columns = renamed_features

		feature_dfs.append(df)

	# check shared samples
	if not shared_labels:
		raise ValueError("No shared labels found across the feature datasets.")
	shared_labels = sorted(shared_labels)
	feature_dfs = [df.loc[shared_labels] for df in feature_dfs]

	# combine features
	X = pd.concat(feature_dfs, axis=1).to_numpy().astype(float)

	# remove columns with NaN
	if remove_nan:
		cols_with_nan = np.isnan(X).any(axis=0)
		X = X[:, ~cols_with_nan]
		feature_names = [name for (name, keep) in zip(feature_names, ~cols_with_nan) if keep]

	if verbose:
		print(f' - After combining types and removing NaN: {X.shape}')

	# remove highly correlated columns
	if X.shape[1] > 1:
		correlation_matrix = np.corrcoef(X, rowvar=False)
		correlated_columns = set()
		for i in range(len(correlation_matrix)):
			for j in range(i + 1, len(correlation_matrix)):
				if abs(correlation_matrix[i, j]) > correlation_threshold:
					correlated_columns.add(j)
		uncorrelated_columns = [i for i in range(X.shape[1]) if i not in correlated_columns]
		X = X[:, uncorrelated_columns]
		feature_names = [feature_names[i] for i in uncorrelated_columns]
		if verbose:
			print(f' - After removing columns with correlation > {correlation_threshold}: {X.shape}')

	# Process Y and drop rows with missing targets
	if responses is not None:
		Y = pd.concat(Y_df_list, axis=1).loc[shared_labels].to_numpy()
		rows_with_nan = np.isnan(Y).any(axis=1)
		X = X[~rows_with_nan, :]
		Y = Y[~rows_with_nan, :]
	else:
		Y = None

	if verbose:
		for category in feature_types:
			surviving = [name for name in feature_names if name in category_feature_names[category]]
			print(f" - Final {category} features: {len(surviving)}")
		print(f"\n - Final covariate dimensions: {X.shape}")
		if Y is not None:
			print(f" - Final target dimensions: {Y.shape}")

	return {'X': X, 'Y': Y, 'feature_names': feature_names, 'responses': responses}

#--------------------------------
# Helpers
#--------------------------------
def load_dataframe(file_path):
	df = pd.read_csv(file_path, sep='\t').T
	df.columns = df.iloc[0]
	df = df.drop(df.index[0])
	df.index.name = 'attrib_name'
	return df

def process_responses(responses, file_paths, verbose=False):
	shared_labels = None
	Y_df_list = []
	for response_type, response_name in responses:
		response_df = load_dataframe(file_paths[response_type])
		replace_map = {
			'gender': {'male': 0, 'female': 1},
			'histological_type': {'infiltratingductalcarcinoma': 0, 'infiltratinglobularcarcinoma': 1},
			'pathologic_stage': {'stagei': 1, 'stageii': 2, 'stageiii': 3, 'stageiv': 4},
			'pathology_M_stage': {'m0': 0, 'm1': 1},
			'pathology_N_stage': {'n0': 0, 'n1': 1, 'n2': 2, 'n3': 3},
			'pathology_T_stage': {'t1': 1, 't2': 2, 't3': 3, 't4': 4},
			'radiation_therapy': {'no': 0, 'yes': 1}
		}
		if response_name in replace_map:
			response_df[response_name] = (
				response_df[response_name]
				.replace(replace_map[response_name])
				.infer_objects(copy=False)
			)
		Y = response_df[response_name].apply(pd.to_numeric, errors='coerce')
		Y_df_list.append(Y)
		if shared_labels is None:
			shared_labels = set(Y.index)
		else:
			shared_labels = shared_labels.intersection(Y.index)
	return Y_df_list, shared_labels

def standardize_feature_names(category, df_columns, feature_types):
	if category == 'mirna':
		renamed = [
			name.replace('hsa-mir-', 'miR-').replace('hsa-let-', 'let-')
			for name in df_columns
		]
	elif category == 'rppa':
		renamed = [
			f"{name.split('|')[0]}\n({name.split('|')[1]})" if '|' in name else name
			for name in df_columns
		]
		if 'rnaseq' in feature_types:
			renamed = [f"{name}#" for name in renamed]
	elif category == 'rnaseq':
		renamed = df_columns
		if 'rppa' in feature_types:
			renamed = [f"{name}*" for name in renamed]
	else:
		renamed = df_columns
	return renamed




