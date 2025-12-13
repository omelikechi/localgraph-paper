# Generate graphical model data

import numpy as np
from sklearn.preprocessing import StandardScaler

import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

#--------------------------------
# Generate data
#--------------------------------
def block_graph(n, lmin, lmax, block_sizes, block_degree, block_magnitude, 
		connector_degree, connector_magnitude, sigma=0, random_seed=None):

	if random_seed is not None:
		np.random.seed(random_seed)

	n_blocks = len(block_sizes)

	Omega = generate_block_theta(block_sizes, block_degree, block_magnitude, connector_degree, connector_magnitude, sigma)
	Omega = make_posdef(Omega, lmin, lmax)
	p = Omega.shape[0]

	A = np.copy(Omega)
	A[A != 0] = 1
	np.fill_diagonal(A,0)

	Sigma = np.linalg.inv(Omega)
	X = np.random.multivariate_normal(np.zeros(p), Sigma, size=n)
	X = StandardScaler().fit_transform(X)

	# compute maximum correlation with response(s)
	Sigma_hat = X.T @ X / n
	abs_cor = np.abs(Sigma_hat)
	np.fill_diagonal(abs_cor,0)
	max_cor_response = []
	for i in range(block_sizes[0]):
		max_cor = np.max(abs_cor[i,:])
		max_cor_response.append(max_cor)

	target_features = np.arange(block_sizes[0])

	return {'X':X, 'A':A.astype(int), 'feature_names':None, 'target_features':target_features, 'max_cor_response':max_cor_response}

#--------------------------------
# Helpers
#--------------------------------
# generate precision matrix for the block design
def generate_block_theta(block_sizes, block_degree, block_magnitude, connector_degree, connector_magnitude, sigma):
	n_blocks = len(block_sizes)
	p = sum(block_sizes)
	block_sparsity = [block_degree[i] / block_sizes[i] for i in range(n_blocks)]
	connector_sparsity = [connector_degree[i] / max(block_sizes[i],block_sizes[i+1]) for i in range(n_blocks-1)]
	Omega = np.zeros((p, p))
	start_idx = 0
	# handle the block matrices
	for i in range(n_blocks):
		end_idx = start_idx + block_sizes[i]
		# fill block matrix
		for row in range(start_idx, end_idx):
			for col in range(row, end_idx):
				if np.random.rand() <= block_sparsity[i]:
					sign = np.random.choice([-1, 1])
					Omega[row, col] = sign * np.random.normal(block_magnitude[i], sigma)
		if i < n_blocks - 1:
			# handle connector matrices between block i and block i+1
			next_block_start = end_idx
			next_block_end = next_block_start + block_sizes[i + 1]
			for row in range(start_idx, end_idx):  # Rows of block i
				for col in range(next_block_start, next_block_end):  # Columns of block i+1
					if np.random.rand() <= connector_sparsity[i]:
						sign = np.random.choice([-1, 1])
						Omega[row, col] = sign * np.random.normal(connector_magnitude[i], sigma)
		# move start index to the next block start
		start_idx = end_idx
	# symmetrize matrix
	Omega = Omega + Omega.T - np.diag(np.diag(Omega))
	return Omega

# ensure precision matrix is positive definite
def make_posdef(Omega, lmin, lmax):
	p = Omega.shape[0]

	eigenvalues = np.linalg.eigvalsh(Omega)
	lambda_min = np.min(eigenvalues)
	Omega -= lambda_min * np.eye(p)

	eigenvalues = np.linalg.eigvalsh(Omega)
	lambda_max = np.max(eigenvalues)
	Omega *= ((lmax - lmin) / lambda_max)
	Omega += lmin * np.eye(p)
	
	return Omega


