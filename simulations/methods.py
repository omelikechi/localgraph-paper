# Methods for graphical inference

# from dagma.linear import DagmaLinear
# from dagma.nonlinear import DagmaMLP, DagmaNonlinear
import numpy as np

import rpy2.robjects as robjects
from rpy2.robjects import FloatVector, numpy2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

# import R packages
huge = importr('huge')
mgm = importr('mgm')
silggm = importr('SILGGM')

# #--------------------------------
# # dagma
# #--------------------------------
# def dagma(X, linear=True, lambda1=0.02, lambda2=0.005):
# 	"""
# 	Estimate the DAG using the DAGMA method (linear or nonlinear).
	
# 	Parameters:
# 	- X: Input data (numpy array).
# 	- method: 'linear' or 'nonlinear' to select the model type.
# 	- lambda1: L1 regularization coefficient.
# 	- lambda2: L2 regularization coefficient (only for nonlinear model).
	
# 	Returns:
# 	- Adjacency matrix (binary) estimated by the DAGMA method.
# 	"""
# 	if linear:
# 		model = DagmaLinear(loss_type='l2')  # Linear DAG learning
# 		W_est = model.fit(X, lambda1=lambda1)  # Fit the model
# 	else:
# 		d = X.shape[1]
# 		eq_model = DagmaMLP(dims=[d, 10, 1], bias=True)  # Structural equations
# 		model = DagmaNonlinear(eq_model)  # Nonlinear DAG learning
# 		W_est = model.fit(X, lambda1=lambda1, lambda2=lambda2)  # Fit the model

# 	adjacency_matrix = (W_est != 0).astype(int)  # Convert to binary adjacency matrix
# 	return adjacency_matrix

#--------------------------------
# gipss
#--------------------------------
def run_gipss(X, target_features, fdr_path, fdr_local=None, max_radius=None, feature_names=None, selector='gb'):
	result = gipss(X, target_features, fdr_path, fdr_local=fdr_local, max_radius=max_radius, 
		feature_names=feature_names, selector=selector)
	adjacency_matrix = result['Q_adj']
	return adjacency_matrix

#--------------------------------
# huge
#--------------------------------
def hugeR(X, method, **huge_args):

	with localconverter(robjects.default_converter + numpy2ri.converter):
		X_r = robjects.conversion.py2rpy(X)

	criterion = huge_args.pop('criterion', 'ric')
	lambda_ = huge_args.pop('lambda_', None)

	# apply nonparanormal transformation
	apply_npn = huge_args.pop('apply_npn', False)
	if apply_npn:
		npn = robjects.r['huge.npn']
		X_r = npn(X_r, verbose=False)

	if lambda_ is not None:
		lambda_ = FloatVector([float(lambda_)])
		result = huge.huge(X_r, method=method, verbose=False, **{'lambda': lambda_}, **huge_args)
		A = result.rx2('path')[0]
	else:
		result = huge.huge(X_r, method=method, verbose=False, **huge_args)
		select = huge.huge_select(result, criterion=criterion, verbose=False)
		A = select.rx2('refit')

	with localconverter(robjects.default_converter + numpy2ri.converter):
		A_numpy = np.array(robjects.r['as.matrix'](A))

	adjacency_matrix = (A_numpy != 0).astype(int)

	return adjacency_matrix

#--------------------------------
# mgm
#--------------------------------
def mgmR(X, feature_type=None, level=None):
	if feature_type is None:
		feature_type = ['g'] * X.shape[1]  # Default to Gaussian variables
	if level is None:
		level = [1] * X.shape[1]  # Default level for Gaussian variables

	with localconverter(robjects.default_converter + numpy2ri.converter):
		r_matrix = X

		# Call the R mgm function
		mgm = robjects.r['mgm']
		robjects.r('sink("/dev/null")')
		try:
			result = mgm(data=r_matrix, type=feature_type, level=level)
		finally:
			robjects.r('sink()')

		# Extract the weighted adjacency matrix (equivalent to result$pairwise$wadj in R)
		A = result.rx2('pairwise').rx2('wadj')

		# Convert to NumPy array
		with localconverter(robjects.default_converter + numpy2ri.converter):
			A_numpy = np.array(A)

		# Convert to binary adjacency matrix (1 if nonzero, else 0)
		adjacency_matrix = (A_numpy != 0).astype(int)

	return adjacency_matrix

# #--------------------------------
# # pc
# #--------------------------------
# from causallearn.search.ConstraintBased.PC import pc

# def pc_algorithm(X, alpha=0.05):
# 	"""
# 	Estimate the DAG using the PC algorithm.

# 	Parameters:
# 	- X: Input data (numpy array).
# 	- alpha: Significance level for conditional independence tests (default is 0.05).

# 	Returns:
# 	- Adjacency matrix (binary) estimated by the PC algorithm.
# 	"""
# 	# Run the PC algorithm
# 	pc_graph = pc(X, alpha=alpha, show_progress=False)

# 	# Extract the adjacency matrix directly from the GeneralGraph object
# 	adjacency_matrix = pc_graph.G.graph

# 	# Ensure binary output (0 or 1)
# 	adjacency_matrix = (adjacency_matrix != 0).astype(int)

# 	return adjacency_matrix

#--------------------------------
# silggm
#--------------------------------
def silggmR(X, method, **silggm_args):
	with localconverter(robjects.default_converter + numpy2ri.converter):
		X_r = robjects.conversion.py2rpy(X)

	# Apply nonparanormal transformation if requested
	apply_npn = silggm_args.pop('apply_npn', False)
	if apply_npn:
		npn = robjects.r['huge.npn']
		X_r = npn(X_r, verbose=False)

	silggm_args.setdefault('alpha', 0.05)
	silggm_args.setdefault('global', True)

	# Ensure alpha is an R numeric vector
	alpha = silggm_args['alpha']
	if isinstance(alpha, (list, tuple, np.ndarray)):
		silggm_args['alpha'] = FloatVector(alpha)
		alphas = list(alpha)
	else:
		alphas = [alpha]

	robjects.r('sink("/dev/null")')
	try:
		result = silggm.SILGGM(X_r, method=method, **silggm_args)
	finally:
		robjects.r('sink()')

	global_decision = result.rx2('global_decision')

	adjacency_matrices = {}
	for alpha, mat in zip(alphas, global_decision):
		A = robjects.r['as.matrix'](mat)
		with localconverter(robjects.default_converter + numpy2ri.converter):
			A_numpy = np.array(A)
		adjacency_matrices[float(alpha)] = (A_numpy != 0).astype(int)

	if len(alphas) == 1:
		return adjacency_matrices[alphas[0]]
	else:
		return adjacency_matrices





