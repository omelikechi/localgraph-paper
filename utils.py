# Utility functions for localgraph-paper repository

import numpy as np

def max_cor_response(X, target_features):
	"""
	Compute the maximum absolute correlation of each target feature with all other features.

	Parameters
	--------------------------------
	X : numpy.ndarray
		Data matrix of shape (n, p), where n is the number of samples and p is the number of features.
	target_features : list of int
		Indices of the target features for which correlations are computed.

	Returns
	--------------------------------
	max_cors : list of float
		List of maximum absolute correlations, one value for each target feature in `target_features`.
	"""
	n = X.shape[0]
	Sigma_hat = X.T @ X / n
	abs_cor = np.abs(Sigma_hat)
	np.fill_diagonal(abs_cor, 0)
	max_cors = []
	for i in target_features:
		max_cor = np.max(abs_cor[i, :])
		max_cors.append(max_cor)
		
	return max_cors


	