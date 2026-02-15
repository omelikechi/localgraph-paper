# Run method from R package 'mgm' (mixed graphical models)
"""
- CRAN: https://cran.r-project.org/web/packages/mgm/mgm.pdf
- Paper: https://www.jstatsoft.org/article/view/v093i08
"""

import time

import numpy as np

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter

# import R package
mgm = importr('mgm')

import numpy as np

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter

# import R package
mgm = importr('mgm')

def run_mgm(X, **mgm_args):
	mgm_args.setdefault('k', 2)
	cat_threshold = mgm_args.pop('cat_threshold', 10)

	n, p = X.shape
	feature_type = []
	level = []

	for j in range(p):
		col = X[:, j]
		unique_vals, counts = np.unique(col, return_counts=True)

		# constant variable
		if len(unique_vals) == 1:
			feature_type.append("g")
			level.append(1)
			continue

		if len(unique_vals) <= cat_threshold:
			# if every category has > 1 observation, treat as categorical
			if np.all(counts >= 2):
				feature_type.append("c")
				level.append(len(unique_vals))
			else:
				# fallback to Gaussian
				feature_type.append("g")
				level.append(1)
		else:
			feature_type.append("g")
			level.append(1)

	with localconverter(robjects.default_converter + numpy2ri.converter):
		X_r = robjects.conversion.py2rpy(X)

	mgm_fun = robjects.r['mgm']
	# robjects.r('sink("/dev/null")')
	try:
		start_time = time.time()
		result = mgm_fun(data=X_r, type=feature_type, level=level, **mgm_args)
		runtime = time.time() - start_time
	finally:
		# robjects.r('sink()')
		pass

	A = result.rx2('pairwise').rx2('wadj')
	with localconverter(robjects.default_converter + numpy2ri.converter):
		A_numpy = np.array(A)

	adjacency_matrix = (A_numpy != 0).astype(int)

	return {'adjacency_matrix':adjacency_matrix, 'runtime':runtime}


