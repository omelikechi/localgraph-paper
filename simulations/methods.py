# Methods for graphical inference

import warnings

# from dagma.linear import DagmaLinear
# from dagma.nonlinear import DagmaMLP, DagmaNonlinear
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter

# enable automatic numpy to rpy2 conversion
numpy2ri.activate()

# import R packages
flare = importr('flare')
huge = importr('huge')
mgm = importr('mgm')
silggm = importr('SILGGM')
wgcna = importr('WGCNA')

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
# flare
#--------------------------------
def flareR(X, method='tiger', lambda_=None, criterion='stars', sym='or'):
	r_matrix = numpy2ri.py2rpy(X)

	if lambda_ is None:
		result = flare.sugm(r_matrix, method=method, sym=sym, verbose=False)
		select = flare.sugm_select(result, criterion=criterion, verbose=False)
		A = select.rx2('refit')
	else:
		result = flare.sugm(r_matrix, method=method, **{'lambda': lambda_}, sym=sym, verbose=False)
		A = result.rx2('path')[0]

	A = robjects.r['as.matrix'](A)
	with localconverter(robjects.default_converter + numpy2ri.converter):
		A_numpy = np.array(A)
	adjacency_matrix = (A_numpy != 0).astype(int)

	return adjacency_matrix

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
def hugeR(X, method='glasso', lambda_=None, criterion='ric', sym='or'):
	r_matrix = numpy2ri.py2rpy(X)

	if lambda_ is None:
		result = huge.huge(r_matrix, method=method, sym=sym, verbose=False)
		select = huge.huge_select(result, criterion=criterion, verbose=False)
		A = select.rx2('refit')
	else:
		result = huge.huge(r_matrix, method=method, **{'lambda': lambda_}, verbose=False)
		A = result.rx2('path')[0]

	A = robjects.r['as.matrix'](A)
	with localconverter(robjects.default_converter + numpy2ri.converter):
		A_numpy = np.array(A)
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

	r_matrix = numpy2ri.py2rpy(X)

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
def silggmR(X, method='GFC_L', lambda_=None, alpha=None):
	r_matrix = numpy2ri.py2rpy(X)

	if lambda_ is None:
		lambda_ = robjects.NULL
	if alpha is None:
		alpha = 0.1

	robjects.r('sink("/dev/null")')
	try:
		result = silggm.SILGGM(r_matrix, method=method, alpha=alpha, **{'lambda': lambda_}, **{'global': True})
	finally:
		robjects.r('sink()')

	A = result.rx2('global_decision')
	A = robjects.r['as.matrix'](A)
	with localconverter(robjects.default_converter + numpy2ri.converter):
		A_numpy = np.array(A)
	adjacency_matrix = (A_numpy != 0).astype(int)

	return adjacency_matrix[0]

#--------------------------------
# wgcna
#--------------------------------
def wgcnaR(X, power=6, type="unsigned", threshold=None):
	r_matrix = numpy2ri.py2rpy(X)

	if threshold is None:
		threshold = 0.5

	A = wgcna.adjacency(r_matrix, power=power, type=type)
	with localconverter(robjects.default_converter + numpy2ri.converter):
		A_numpy = np.array(A)
	adjacency_matrix = (A_numpy >= threshold).astype(int)
	for i in range(adjacency_matrix.shape[0]):
		adjacency_matrix[i,i] = 0

	return adjacency_matrix





