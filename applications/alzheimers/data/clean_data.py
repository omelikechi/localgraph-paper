# Clean and save Alzheimers disease single cell data
"""
- Downloaded from https://adsn.ddnetbio.com
- Zenodo: https://zenodo.org/records/17302976
- Paper: "A single-cell atlas of entorhinal cortex from individuals with Alzheimer’s 
	disease reveals cell-type-specific gene expression regulation" by Grubman et al
- DOI: doi:10.1038/s41593-019-0539-4
""" 

import numpy as np
import pandas as pd
import pickle
from scipy.io import mmread

################################
save_data = False
################################

print('Loading raw data...')

X = mmread('adsn_matrix.mtx').tocsr()   # genes x cells
genes = pd.read_csv('adsn_features.tsv', header=None, sep='\t')[0].values
cells = pd.read_csv('adsn_barcodes.tsv', header=None)[0].values

X = X.T

meta = pd.read_csv('adsn_metadata.txt', sep='\t')
meta.index = meta.iloc[:,0]
meta = meta.loc[cells]

print('\nAD vs Control counts:')
print(meta['batchCond'].value_counts())

print('\nCell type counts:')
print(meta['cellType'].value_counts())

#----------------------------------------------------------------
# Gene filtering
#----------------------------------------------------------------
expression_cutoff = 50
gene_detect = np.array((X > 0).sum(axis=0)).ravel()
keep = gene_detect >= expression_cutoff
X = X[:, keep]
genes = genes[keep]

print('After filtering:', X.shape)

#----------------------------------------------------------------
# Log-normalize
#----------------------------------------------------------------
lib_sizes = np.array(X.sum(axis=1)).ravel()
X = X.multiply(1e4 / lib_sizes[:, None])
X.data = np.log1p(X.data)
X = X.tocsr()

#----------------------------------------------------------------
# Save cleaned dataset
#----------------------------------------------------------------
if save_data:
	print('Saving cleaned dataset...')
	out = {'X': X, 'genes': genes, 'meta': meta}
	with open('cleaned_data.pkl', 'wb') as f:
		pickle.dump(out, f)
	print('Saved to cleaned_data.pkl')


