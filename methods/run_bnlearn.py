# Run methods from the R package 'bnlearn'
"""
There are two main functions: run_bnlearn and run_bnlearn_local. The former estimates
the full graph; the latter estimates the local graph of some user-specified radius around
certain target features. In bnlearn, there are two functions for learning the radius-1
neighborhood around a node/feature: learn.mb learns the markov blanket of a node and learn.nbr
learns the parents and children of a node. The only methods compatible with learn.mb are
fast.iamb, gs, iamb, iamb.fdr, and inter.iamb. The only methods compatible with learn.nbr
are hpc, mmpc, pc.stable, and si.hiton.pc.
"""

import time
import numpy as np

import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

bnlearn = importr('bnlearn')

constraint_based = {
	'gs',
	'iamb',
	'fast.iamb',
	'inter.iamb',
	'iamb.fdr',
	'pc.stable',

	# undirected only
	'hpc',
	'mmpc',
	'si.hiton.pc'
}

pairwise_mi = {
	'aracne'
}

def _sanitize_kwargs(d):
	out = {}
	for k, v in d.items():
		if isinstance(v, np.generic):
			out[k] = v.item()
		else:
			out[k] = v
	return out

#----------------------------------------------------------------
# Global version of bnlearn methods
#----------------------------------------------------------------
def run_bnlearn(X, method, **bnlearn_args):
	n, p = X.shape
	adjacency = np.zeros((p, p), dtype=int)

	bnlearn_args = _sanitize_kwargs(bnlearn_args)
	if method in constraint_based:
		bnlearn_args.setdefault('undirected', True)

	with localconverter(robjects.default_converter + numpy2ri.converter):
		X_r = robjects.conversion.py2rpy(X)

	X_r = robjects.r['as.data.frame'](X_r)

	# bnlearn requires column names
	robjects.r.colnames(X_r).ro = [f'V{i + 1}' for i in range(p)]

	start_time = time.time()

	py_method = method.replace('.', '_')
	f = getattr(bnlearn, py_method)

	res = f(X_r, robjects.NULL, **bnlearn_args)

	for i in range(p):
		mb = list(res.rx2('nodes')[i].rx2('mb'))
		for v in mb:
			j = int(v[1:]) - 1
			adjacency[i,j] = 1
			adjacency[j,i] = 1

	runtime = time.time() - start_time

	return {'adjacency_matrix':adjacency, 'runtime':runtime}

#----------------------------------------------------------------
# Local version of bnlearn methods
#----------------------------------------------------------------
def run_bnlearn_local(X, method, target_features, radius=1, criterion=None, verbose=False, **bnlearn_args):
	n, p = X.shape
	adjacency = np.zeros((p, p), dtype=int)

	bnlearn_args = _sanitize_kwargs(bnlearn_args)

	with localconverter(robjects.default_converter + numpy2ri.converter):
		X_r = robjects.conversion.py2rpy(X)
	X_r = robjects.r['as.data.frame'](X_r)
	robjects.r.colnames(X_r).ro = [f'V{i + 1}' for i in range(p)]

	start_time = time.time()

	learn_mb  = robjects.r['learn.mb']
	learn_nbr = robjects.r['learn.nbr']

	mb_methods  = {'fast.iamb', 'gs', 'iamb', 'inter.iamb', 'iamb.fdr'}
	nbr_methods = {'hpc', 'mmpc', 'pc.stable', 'si.hiton.pc'}

	if method in mb_methods:
		local_fun = learn_mb
	elif method in nbr_methods:
		local_fun = learn_nbr
	else:
		raise ValueError(f"Method '{method}' is not supported.")

	visited = set()
	if isinstance(target_features, (list, tuple, np.ndarray)):
		frontier = set(int(i) for i in np.asarray(target_features).ravel())
	else:
		frontier = {int(target_features)}

	for r in range(radius):

		if verbose:
			node_iteration = 1
			print(f'current features: {frontier} (radius = {r + 1}/{radius})')

		new_frontier = set()

		for i in frontier:

			if verbose:
				start = time.time()

			if i in visited:
				continue
			visited.add(i)

			node_name = f'V{i+1}'
			neigh = list(local_fun(X_r, node_name, method=method, **bnlearn_args))

			for v in neigh:
				j = int(str(v)[1:]) - 1

				if criterion == 'forward':
					# only accept edges from current layer to unseen nodes
					if j not in visited:
						adjacency[i, j] = adjacency[j, i] = 1
						new_frontier.add(j)

				else:  # union-style (original behavior)
					adjacency[i, j] = adjacency[j, i] = 1
					if j not in visited:
						new_frontier.add(j)

			if verbose:
				runtime = time.time() - start
				print(f' - iteration {node_iteration}/{len(frontier)} ({runtime:.2f} seconds)')
				node_iteration += 1

		frontier = new_frontier

	runtime = time.time() - start_time
	return {'adjacency_matrix': adjacency, 'runtime': runtime}


