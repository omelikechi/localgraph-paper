# Save cleaned county-level data (Environmental/socioeconomic cancer study)

import numpy as np
import pandas as pd
import os

from load_and_clean import load_and_clean

random_seed = 4161932
np.random.seed(random_seed)

#--------------------------
# Parameters
#--------------------------
responses = ['Mortality', 'Incidence']
features_to_drop = ['Pop2023', 'pct_no_eng', 'Longitude', 'Latitude']
states_to_remove = ['AK', 'HI', 'KS', 'IN']
constant_state_threshold = 5
n_redundant = 1500

data_dir = './raw_data'

#--------------------------
# Load and clean data
#--------------------------
data = load_and_clean(
	data_dir=data_dir,
	data_name='eqi2000.csv',
	metadata_file='eqi2000_feature_names.xlsx',
	responses=responses,
	features_to_drop=features_to_drop,
	constant_state_threshold=constant_state_threshold,
	n_redundant=n_redundant,
	remove_nan=True,
	states_to_remove=states_to_remove,
	verbose=True
)

X = data['X']
Y = data['Y']
feature_names = data['feature_names']

# Remove Union County, Florida (extreme outlier in both targets)
idx = responses.index('Incidence')  # fix: index from responses, not feature_names
mask = Y[:, idx] <= 1000
X = X[mask]
Y = Y[mask]

print(f'\nFinal shape after removing Union County, FL: X = {X.shape}, Y = {Y.shape}')

#--------------------------
# Save cleaned full data
#--------------------------
out_dir = './cleaned_data'
os.makedirs(out_dir, exist_ok=True)

full_data = np.column_stack((Y, X))
full_columns = feature_names

df_full = pd.DataFrame(full_data, columns=full_columns)
df_full.to_csv(os.path.join(out_dir, 'cleaned_data.csv'), index=False)
