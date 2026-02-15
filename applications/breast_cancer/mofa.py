# Apply MOFA+ to breast cancer data

import matplotlib.pyplot as plt
from mofapy2.run.entry_point import entry_point
import numpy as np
import pandas as pd

# import data
data_path = 'data/cleaned_data/cleaned_data.csv'
df = pd.read_csv(data_path)

target_names = ['histological_type', 'pathologic_stage', 'status']

# split into distinct modalities (views)
rna_cols = []
mirna_cols = []
rppa_cols = []
other_cols = []

for col in df.columns:
	if col in target_names:
		continue
	if 'miR-' in col or 'let-' in col:
		mirna_cols.append(col)
	elif '#' in col:
		rppa_cols.append(col)
	elif '*' in col:
		rna_cols.append(col)

# clinical encoding
binary_targets = ['histological_type', 'status']
stage_dummies = pd.get_dummies(df['pathologic_stage'], prefix='stage')

clinical_df = pd.concat([df[binary_targets], stage_dummies], axis=1)
clinical_cols = clinical_df.columns.tolist()

# views
views = {
	'RNA': [df[rna_cols].values],
	'miRNA': [df[mirna_cols].values],
	'RPPA': [df[rppa_cols].values],
	'Clinical': [clinical_df.values]
}
likelihoods = ['gaussian', 'gaussian', 'gaussian', 'bernoulli']

# prin dimensions
print('MOFA input dimensions:')
for k, v in views.items():
	Xk = v[0]
	print(f'  {k:10s}: n = {Xk.shape[0]}, p = {Xk.shape[1]}')

# run MOFA+
ent = entry_point()
ent.set_data_matrix(data=views, likelihoods=likelihoods)
ent.set_model_options(factors=10)
ent.set_train_options(seed=1)
ent.build()
ent.run()
model = ent.model

# extract clinical loadings
W = model.nodes['W'].getExpectation()
W_clinical = W[3]

print('Clinical W shape:', W_clinical.shape)

summary = []

# binary targets
for i, name in enumerate(['histological_type', 'status']):
	loadings = W_clinical[i,:]
	k = np.argmax(np.abs(loadings))
	summary.append({
		'target': name,
		'best_factor': f'Factor{k+1}',
		'abs_loading': float(np.abs(loadings[k]))
	})

# pathologic stage (one-hot encoding)
stage_idx = list(range(2,6))
stage_loadings = W_clinical[stage_idx, :].mean(axis=0)
k = np.argmax(np.abs(stage_loadings))

summary.append({
	'target': 'pathologic_stage',
	'best_factor': f'Factor{k+1}',
	'abs_loading': float(np.abs(stage_loadings[k]))
})

print(pd.DataFrame(summary))

# Visualize clinical loadings
factor_names = [f'Factor{k+1}' for k in range(W_clinical.shape[1])]

fig, ax = plt.subplots(figsize=(7,4))

# Histological type
ax.plot(np.abs(W_clinical[0, :]), 'o-', label='Histological type')

# Status
ax.plot(np.abs(W_clinical[1, :]), 'o-', label='Status')

# Pathologic stage (average across one-hot columns)
ax.plot(np.abs(stage_loadings), 'o-', label='Pathologic stage')

ax.set_xticks(range(len(factor_names)))
ax.set_xticklabels(factor_names, rotation=45)
ax.set_ylabel('Absolute loading')
ax.set_title('MOFA clinical-variable loadings across latent factors')
ax.legend(frameon=False)

fig.tight_layout()
plt.show()


