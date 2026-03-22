# Collection of methods for computing q-values

import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter


# q-values from knockoffs (R package)
def knockoff_qvalues(X, y, alpha_list, stat='glmnet_coefdiff', mu=None, Sigma=None):
	knockoff_r_code = f"""
	suppressMessages({{library(knockoff)}})
	run_knockoff_filter <- function(X, y, alpha_list, mu, Sigma) {{
		if (!is.null(mu) && !is.null(Sigma)) {{
			knockoffs = function(X) create.gaussian(X, mu, Sigma)
		}} else {{
			knockoffs = function(X) create.second_order(X)
		}}
		k_stat = function(X, Xk, y) stat.{stat}(X, Xk, y, nfolds=5)
		result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat, fdr=alpha_list[1])
		W = result$statistic
		qvals <- rep(1, length(W))
		for (alpha in alpha_list) {{
			t = knockoff.threshold(W, fdr=alpha)
			selected = which(W >= t)
			qvals[selected] = pmin(qvals[selected], alpha)
		}}
		return(qvals)
	}}
	"""
	X = np.asarray(X, dtype=float)
	y = np.asarray(y, dtype=float)
	with localconverter(robjects.default_converter + numpy2ri.converter):
		X_r = robjects.conversion.py2rpy(X)
		y_r = robjects.conversion.py2rpy(y)
		mu_r = robjects.NULL if mu is None else robjects.conversion.py2rpy(mu)
		Sigma_r = robjects.NULL if Sigma is None else robjects.conversion.py2rpy(Sigma)
	robjects.r(knockoff_r_code)
	run_knockoff_filter_r = robjects.globalenv['run_knockoff_filter']
	alpha_list_r = robjects.FloatVector(alpha_list)
	qvals_r = run_knockoff_filter_r(X_r, y_r, alpha_list_r, mu_r, Sigma_r)
	with localconverter(robjects.default_converter + numpy2ri.converter):
		qvals = np.array(qvals_r)
	p = len(qvals)
	q_values = {j: qvals[j] for j in range(p)}

	return {'q_values': q_values}


# Benjamini-Hochberg and Benjamini-Yekutieli
def padjust(X, y, method='bh'):
	method = 'fdr_' + method
	n, p = X.shape
	if n > p:
		X_full = sm.add_constant(X)
		results = sm.OLS(y, X_full).fit()
		p_values = results.pvalues[1:]
	else:
		p_values = np.zeros(p)
		for j in range(p):
			Xi = sm.add_constant(X[:, j])
			results = sm.OLS(y, Xi).fit()
			p_values[j] = results.pvalues[1]
	_, qvals, _, _ = multipletests(p_values, method=method)
	q_values = {j: qvals[j] for j in range(p)}

	return {'q_values': q_values}




























# import numpy as np
# import statsmodels.api as sm
# from statsmodels.stats.multitest import multipletests

# import rpy2.robjects as robjects
# from rpy2.robjects import numpy2ri
# from rpy2.robjects.conversion import localconverter


# # q-values from knockoffs
# def knockoff_qvalues(X, y, alpha_list, stat='glmnet_coefdiff', mu=None, Sigma=None):

# 	knockoff_r_code = f"""
# 	suppressMessages({{library(knockoff)}})

# 	run_knockoff_filter <- function(X, y, alpha_list, mu, Sigma) {{
# 		if (!is.null(mu) && !is.null(Sigma)) {{
# 			knockoffs = function(X) create.gaussian(X, mu, Sigma)
# 		}} else {{
# 			knockoffs = function(X) create.second_order(X)
# 		}}

# 		k_stat = function(X, Xk, y) stat.{stat}(X, Xk, y, nfolds=5)

# 		result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat, fdr=alpha_list[1])
# 		W = result$statistic

# 		qvals <- rep(1, length(W))

# 		for (alpha in alpha_list) {{
# 			t = knockoff.threshold(W, fdr=alpha)
# 			selected = which(W >= t)
# 			qvals[selected] = pmin(qvals[selected], alpha)
# 		}}

# 		return(qvals)
# 	}}
# 	"""

# 	X = np.asarray(X, dtype=float)
# 	y = np.asarray(y, dtype=float)

# 	with localconverter(robjects.default_converter + numpy2ri.converter):
# 		X_r = robjects.conversion.py2rpy(X)
# 		y_r = robjects.conversion.py2rpy(y)
# 		mu_r = robjects.NULL if mu is None else robjects.conversion.py2rpy(mu)
# 		Sigma_r = robjects.NULL if Sigma is None else robjects.conversion.py2rpy(Sigma)

# 	robjects.r(knockoff_r_code)
# 	run_knockoff_filter_r = robjects.globalenv['run_knockoff_filter']

# 	alpha_list_r = robjects.FloatVector(alpha_list)

# 	qvals_r = run_knockoff_filter_r(X_r, y_r, alpha_list_r, mu_r, Sigma_r)

# 	with localconverter(robjects.default_converter + numpy2ri.converter):
# 		qvals = np.array(qvals_r)

# 	p = len(qvals)
# 	q_values = {j: qvals[j] for j in range(p)}

# 	return {'q_values': q_values}

# # Benjamini-Hochberg and Benjamini-Yekutieli
# def padjust(X, y, method='bh'):

# 	method = 'fdr_' + method

# 	n, p = X.shape

# 	if n > p:
# 		X_full = sm.add_constant(X)
# 		results = sm.OLS(y, X_full).fit()
# 		p_values = results.pvalues[1:]
# 	else:
# 		p_values = np.zeros(p)
# 		for j in range(p):
# 			Xi = sm.add_constant(X[:, j])
# 			results = sm.OLS(y, Xi).fit()
# 			p_values[j] = results.pvalues[1]

# 	_, qvals, _, _ = multipletests(p_values, method=method)

# 	q_values = {j: qvals[j] for j in range(p)}

# 	return {'q_values': q_values}


# 	