# Metadata for different graph estimation methods

# method names for plotting, tables, etc
method_names = {
	# bnlearn
	'aracne':'ARACNE',
	'fast_iamb':'FastIAMB',
	'hpc':'HPC',
	'iamb':'IAMB',
	'mmpc':'MMPC',
	'pc_stable':'StablePC',
	'si_hiton_pc':'SIHPC',

	# bnlearn local
	'fast_iamb_local':'FastIAMB(L)',
	'hpc_local':'HPC(L)',
	'iamb_local':'IAMB(L)',
	'mmpc_local':'MMPC(L)',
	'pc_stable_local':'StablePC(L)',
	'si_hiton_pc_local':'SIHPC(L)',

	# huge
	'glasso':'Glasso',
	'mb':'NLasso',

	# pfs
	'pfs':'PFS',

	# silggm
	'bnwsl':'BNWSL',
	'dsgl':'DSGL',
	'dsnwsl':'DSNWSL',
	'gfcl':'GFCL',
	'gfcsl':'GFCSL'
}

bnlearn_methods = ['aracne', 'fast_iamb', 'hpc', 'iamb', 'mmpc', 'pc_stable', 'si_hiton_pc']
bnlearn_local_methods = [method + '_local' for method in bnlearn_methods]
huge_methods = ['glasso', 'mb']
silggm_methods = ['bnwsl', 'dsgl', 'dsnwsl', 'gfcl', 'gfcsl']

def method_type(method):
	if method in bnlearn_methods:
		method_type = 'bnlearn'
	elif method in bnlearn_local_methods:
		method_type = 'bnlearn_local'
	elif method in huge_methods:
		method_type = 'huge'
	elif method in silggm_methods:
		method_type = 'silggm'
	elif 'pfs' in method:
		method_type = method
	else:
		method_type = None
		print(f'Error: Method name {method} not supported.')

	return method_type


