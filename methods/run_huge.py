# Run different methods from the R package 'huge'
"""
- CRAN: https://cran.r-project.org/web/packages/huge/huge.pdf
- Paper: https://jmlr.csail.mit.edu/papers/volume13/zhao12a/zhao12a.pdf
"""

import time

import numpy as np

import rpy2.robjects as robjects
from rpy2.robjects import FloatVector, numpy2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

# import R package
huge = importr('huge')

def run_huge(X, method, **huge_args):

	with localconverter(robjects.default_converter + numpy2ri.converter):
		X_r = robjects.conversion.py2rpy(X)

	criterion = huge_args.pop('criterion', 'ric')
	lambda_ = huge_args.pop('lambda_', None)

	# Apply nonparanormal transformation if requested
	apply_npn = huge_args.pop('apply_npn', False)
	if apply_npn:
		npn = robjects.r['huge.npn']
		X_r = npn(X_r, verbose=False)

	if lambda_ is not None:
		lambda_ = FloatVector([float(lambda_)])
		start_time = time.time()
		result = huge.huge(X_r, method=method, verbose=False, **{'lambda': lambda_}, **huge_args)
		runtime = time.time() - start_time
		A = result.rx2('path')[0]
	else:
		start_time = time.time()
		result = huge.huge(X_r, method=method, verbose=False, **huge_args)
		select = huge.huge_select(result, criterion=criterion, verbose=False)
		runtime = time.time() - start_time
		A = select.rx2('refit')

	with localconverter(robjects.default_converter + numpy2ri.converter):
		A_numpy = np.array(robjects.r['as.matrix'](A))

	adjacency_matrix = (A_numpy != 0).astype(int)

	return {'adjacency_matrix':adjacency_matrix, 'runtime':runtime}


