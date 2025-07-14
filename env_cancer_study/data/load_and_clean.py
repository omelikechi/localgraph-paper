# Functions for loading county-level environmental, socioeconomic, demographic, and health data

import os
import numpy as np
import pandas as pd

def load_and_clean(
	data_dir,
	data_name='eqi2000.csv',
	metadata_file='eqi2000_feature_names.xlsx',
	responses=None,
	features_to_drop=None,
	constant_state_threshold=None,
	n_redundant=None,
	remove_nan=True,
	states_to_remove=None,
	verbose=False
):
	if responses is None:
		raise ValueError("Please specify the response variables.")
	features_to_drop = features_to_drop or []
	states_to_remove = states_to_remove or []

	# Load raw data
	raw_df = pd.read_csv(os.path.join(data_dir, data_name))
	if verbose:
		print(f'Original dataset size: {raw_df.shape}')

	# Identify and drop features constant within many states
	if constant_state_threshold is not None:
		constant_feats = identify_constant_features(raw_df, states_to_remove=states_to_remove)
		for feat, (states, count) in constant_feats.items():
			if count >= constant_state_threshold:
				if verbose:
					print(f"Feature '{feat}' is constant in {count} states: {states}")
				features_to_drop.append(feat)

	# Preprocess
	df = preprocess_eqi_df(raw_df, responses, states_to_remove)

	# Drop redundant features with very low variation
	if n_redundant is not None:
		redundant = [
			col for col in df.drop(columns=responses)
			if df[col].value_counts(dropna=False).max() > n_redundant
		]
		if verbose and redundant:
			print(f" - {len(redundant)} redundant features dropped (>{n_redundant} identical values).")
		features_to_drop.extend(redundant)

	# Drop explicitly listed features
	df.drop(columns=[col for col in features_to_drop if col in df.columns], inplace=True, errors='ignore')

	# Extract response and covariate matrices
	Y = df[responses].to_numpy()
	X = df.drop(columns=responses).to_numpy(dtype=float)

	# Remove rows with NaN if requested
	if remove_nan:
		valid = ~np.isnan(X).any(axis=1)
		X, Y = X[valid], Y[valid]

	# Map human-readable feature names
	feature_names = responses + df.drop(columns=responses).columns.tolist()
	mapping_df = pd.read_excel(os.path.join(data_dir, metadata_file), engine='openpyxl', dtype=str)
	mapping_dict = dict(zip(mapping_df["Variable Name"].str.strip(), mapping_df["Updated Variable Name"].str.strip()))
	feature_names = [mapping_dict.get(name, name) for name in feature_names]

	if verbose:
		print(f'Final dataset dimensions: X={X.shape}, Y={Y.shape}')

	return {
		'X': X,
		'Y': Y,
		'feature_names': feature_names,
		'responses': responses,
		'dropped_features': features_to_drop,
		'index': df.index
	}


def preprocess_eqi_df(df, responses, states_to_remove):
	df = df.copy()

	# Drop unwanted states
	if states_to_remove:
		df = df[~df["State"].isin(states_to_remove)]

	# Parse and clean FIPS
	if 'FIPS' in df.columns:
		df['FIPS'] = pd.to_numeric(df['FIPS'].astype(str).str.extract(r'(\d+)')[0], errors='coerce')
		df.dropna(subset=['FIPS'], inplace=True)
		df['FIPS'] = df['FIPS'].astype(int)

	# Standardize column names
	df.columns = df.columns.str.strip().str.replace('\xa0', ' ', regex=True).str.replace('\u200b', '', regex=True)

	# Drop unnecessary metadata columns
	df.drop(columns=['State', 'County_Name', 'County'], errors='ignore', inplace=True)

	# Convert all to numeric
	df = df.apply(pd.to_numeric, errors='coerce')

	# Drop rows missing any response
	df.dropna(subset=responses, inplace=True)

	# Drop rows that are all NaN
	df.dropna(how='all', inplace=True)

	# Add jitter to certain categorical-like features
	for noisy_col in ['RadonZone', 'cat_RUCC']:
		if noisy_col in df.columns:
			df[noisy_col] += np.random.normal(0, 0.05, size=df.shape[0])

	df.set_index('FIPS', inplace=True)
	return df

def identify_constant_features(df, state_col='State', states_to_remove=None):
	constant_features = {}
	feature_cols = [col for col in df.columns if col not in ['FIPS', state_col, 'County_Name']]
	for feat in feature_cols:
		constant_states = []
		for state, group in df.groupby(state_col):
			if state == 'DC' or (states_to_remove and state in states_to_remove):
				continue
			if group[feat].nunique() == 1:
				constant_states.append(state)
		if constant_states:
			constant_features[feat] = (constant_states, len(constant_states))
	return constant_features
