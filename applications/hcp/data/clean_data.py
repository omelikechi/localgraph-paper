# Prepare cleaned Human Connectome Project (HCP) Young Adult dataset for analysis
"""
Raw data available here: https://balsa.wustl.edu/project?project=HCP_YA

This script constructs the final analysis-ready dataset used in the paper.

Overview
--------------------------------
- target: NIH Toolbox Fluid Cognition Composite score
- covariates: whole-brain structural morphometry (FreeSurfer) + summary-level phenotypes
- input: raw HCP Young Adult subject-level CSV + official HCP data dictionary
- output: cleaned numeric design matrix with documented variable selection and imputation
- data can be age-adjusted or unadjusted; the paper uses age-adjusted
"""

import pickle
import numpy as np
import pandas as pd
import re

#----------------------------------------------------------------
# Setup
#----------------------------------------------------------------

save_data = False

# toggle for age-adjusted vs unadjusted phenotypes
age_adjusted = True
add_on = 'AgeAdj' if age_adjusted else 'Unadj'

# load raw hcp data
data_path = 'hcp_ya_subjects.csv'
df = pd.read_csv(data_path)

# load hcp data dictionary
data_dict = pd.read_csv('HCP_S1200_DataDictionary_Oct_30_2023.csv')

# target variable
target = f'CogFluidComp_{add_on}'
target_names = [target]

#----------------------------------------------------------------
# Drop cognition variables except target
#----------------------------------------------------------------
cog_mask = data_dict['category'].str.contains('Cognition', case=False, na=False)
cog_vars = set(data_dict.loc[cog_mask, 'columnHeader'])
drop_cog = [c for c in df.columns if c in cog_vars and c != target]
df = df.drop(columns=drop_cog)
print(f'Dropped {len(drop_cog)} cognition variables')

#----------------------------------------------------------------
# Drop bookkeeping and metadata variables
#----------------------------------------------------------------
explicit_drop = ['Subject', 'Release', 'Acquisition', 'Gender', 'Age', 'QC_Issue']
pattern_drop = re.compile(r'(_Compl$|_PctCompl$|_Count$|_Avail$|ReconVrs|Scanner|Session)', re.IGNORECASE)

df = df.drop(columns=[c for c in df.columns if c in explicit_drop or pattern_drop.search(c)])

#----------------------------------------------------------------
# Enforce age-adjusted or unadjusted variables globally
#----------------------------------------------------------------
if age_adjusted:
	df = df.drop(columns=[c for c in df.columns if c.endswith('_Unadj')])
else:
	df = df.drop(columns=[c for c in df.columns if c.endswith('_AgeAdj')])

#----------------------------------------------------------------
# Retain summary-level phenotypes only
#----------------------------------------------------------------
wm_keep = {'WM_Task_Acc'}
emotion_keep = {'Emotion_Task_Acc'}
language_keep = {'Language_Task_Acc'}

motor_keep = {
	f'Dexterity_{add_on}',
	f'Endurance_{add_on}',
	f'Strength_{add_on}',
	'GaitSpeed_Comp'
}

sensory_keep = {
	f'Odor_{add_on}',
	f'Taste_{add_on}',
	'Mars_Log_Score',
	'Noise_Comp'
}

personality_keep = {'NEOFAC_A', 'NEOFAC_C', 'NEOFAC_E', 'NEOFAC_N', 'NEOFAC_O'}

#----------------------------------------------------------------
# Final column selection
#----------------------------------------------------------------
cols_keep = []
for c in df.columns:
	if c in target_names:
		cols_keep.append(c)
	elif c.startswith('FS_'):
		cols_keep.append(c)
	elif c in wm_keep:
		cols_keep.append(c)
	elif c in emotion_keep:
		cols_keep.append(c)
	elif c in language_keep:
		cols_keep.append(c)
	elif c in motor_keep:
		cols_keep.append(c)
	elif c in sensory_keep:
		cols_keep.append(c)
	elif c in personality_keep:
		cols_keep.append(c)

df = df.loc[:, cols_keep]
print(f'Kept {len(cols_keep)} variables after phenotype consolidation')

#----------------------------------------------------------------
# Clean numeric data matrix
#----------------------------------------------------------------
# retain numeric columns only
df = df.select_dtypes(include=[np.number])

# drop rows with missing target values only
df = df.dropna(subset=target_names)

# drop zero-variance columns (all non-missing values identical)
# note: nunique(dropna=False) fails when a column is all-zero *except* for NaNs
# use variance instead
variances = df.var(skipna=True)
zero_var_cols = variances[variances == 0].index.tolist()
if zero_var_cols:
	df = df.drop(columns=zero_var_cols)
	print(f'Dropped {len(zero_var_cols)} zero-variance variables')

# mean-impute remaining missing values
for c in df.columns:
	if df[c].isna().any():
		df[c] = df[c].fillna(df[c].mean())

# convert to numpy
X = df.to_numpy()
feature_names = df.columns.tolist()

# identify target indices
target_indices = [feature_names.index(t) for t in target_names]

print(f'\nFinal dataset')
print(f'--------------------------------')
print(f'n = {X.shape[0]}')
print(f'p = {X.shape[1]}')
print(feature_names)

#----------------------------------------------------------------
# Save cleaned dataset
#----------------------------------------------------------------
if save_data:
	out = {
		'X': X,
		'feature_names': feature_names,
		'target': target,
		'target_indices': target_indices,
		'age_adjusted': age_adjusted,
		'n': X.shape[0],
		'p': X.shape[1]
	}

	with open('hcp_cleaned_data.pkl', 'wb') as f:
		pickle.dump(out, f)

	print('Saved cleaned dataset to hcp_cleaned_data.pkl')
