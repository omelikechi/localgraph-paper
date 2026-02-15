# Save the cleaned data using the load_and_clean function from load_and_clean.py

import os
import pandas as pd
from load_and_clean import load_and_clean

# Define settings
feature_types = ['rnaseq', 'mirna', 'rppa']
responses = [('clinical', 'status'), ('clinical', 'histological_type'), ('clinical', 'pathologic_stage')]
expression_threshold = 75
variance_threshold = 75

# Run cleaning
data = load_and_clean(
	feature_types=feature_types,
	responses=responses,
	expression_threshold=expression_threshold,
	variance_threshold=variance_threshold,
	verbose=True
)

# Combine targets and covariates into a single DataFrame
targets = pd.DataFrame(data['Y'], columns=[r[1] for r in responses])
covariates = pd.DataFrame(data['X'], columns=data['feature_names'])
df = pd.concat([targets, covariates], axis=1)

print(df.shape)

# Save as a single cleaned dataset
os.makedirs('./cleaned_data', exist_ok=True)
df.to_csv('./cleaned_data/cleaned_data_var75.csv', index=False)


