# Run different methods from the R package 'huge'
"""
- CRAN: https://cran.r-project.org/web/packages/SILGGM/SILGGM.pdf
- Paper: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006369
- Runs the following methods (refs for each method available at CRAN link above):
	- Bivariate nodewise scaled lasso (B_NW_SL)
	- De-sparsified graphical lasso (D-S_GL)
	- De-sparsified nodewise scaled lasso (D-S_NW_SL)
	- GGM estimation with FDR control
		- Lasso (GFC_L)
		- Scaled lasso (GFC_SL)
"""

import time

import numpy as np

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import FloatVector, numpy2ri
from rpy2.robjects.conversion import localconverter

# import R package
huge = importr('huge')
silggm = importr('SILGGM')

def run_silggm(X, method, **silggm_args):
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
		start_time = time.time()
		result = silggm.SILGGM(X_r, method=method, **silggm_args)
		runtime = time.time() - start_time
	finally:
		robjects.r('sink()')

	global_decision = result.rx2('global_decision')

	adjacency_matrix = {}
	for alpha, mat in zip(alphas, global_decision):
		A = robjects.r['as.matrix'](mat)
		with localconverter(robjects.default_converter + numpy2ri.converter):
			A_numpy = np.array(A)
		adjacency_matrix[float(alpha)] = (A_numpy != 0).astype(int)
	if len(alphas) == 1:
		adjacency_matrix = adjacency_matrix[alphas[0]]

	return {'adjacency_matrix':adjacency_matrix, 'runtime':runtime, 'target_fdrs':alphas}


