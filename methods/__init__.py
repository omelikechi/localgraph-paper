# Main function for running different graph estimation methods

from .metadata import *
from .run_bnlearn import run_bnlearn, run_bnlearn_local
from .run_huge import run_huge
from .run_mgm import run_mgm
from .run_silggm import run_silggm

all_methods = {
	# SILGGM
	'bnwsl':  lambda X, **kw: run_silggm(X, method='B_NW_SL', **kw),
	'dsgl':   lambda X, **kw: run_silggm(X, method='D-S_GL', **kw),
	'dsnwsl': lambda X, **kw: run_silggm(X, method='D-S_NW_SL', **kw),
	'gfcl':   lambda X, **kw: run_silggm(X, method='GFC_L', **kw),
	'gfcsl':  lambda X, **kw: run_silggm(X, method='GFC_SL', **kw),

	# huge
	'glasso': lambda X, **kw: run_huge(X, method='glasso', **kw),
	'mb':     lambda X, **kw: run_huge(X, method='mb', **kw),

	# mgm
	'mgm': run_mgm,

	# bnlearn (global)
	'aracne':      lambda X, **kw: run_bnlearn(X, method='aracne', **kw),
	'fast_iamb':   lambda X, **kw: run_bnlearn(X, method='fast.iamb', **kw),
	'gs':          lambda X, **kw: run_bnlearn(X, method='gs', **kw),
	'hpc':         lambda X, **kw: run_bnlearn(X, method='hpc', **kw),
	'iamb':        lambda X, **kw: run_bnlearn(X, method='iamb', **kw),
	'iamb_fdr':    lambda X, **kw: run_bnlearn(X, method='iamb.fdr', **kw),
	'inter_iamb':  lambda X, **kw: run_bnlearn(X, method='inter.iamb', **kw),
	'mmpc':        lambda X, **kw: run_bnlearn(X, method='mmpc', **kw),
	'pc_stable':   lambda X, **kw: run_bnlearn(X, method='pc.stable', **kw),
	'si_hiton_pc': lambda X, **kw: run_bnlearn(X, method='si.hiton.pc', **kw),

	# bnlearn (local)
	'fast_iamb_local':   lambda X, **kw: run_bnlearn_local(X, method='fast.iamb', **kw),
	'gs_local':          lambda X, **kw: run_bnlearn_local(X, method='gs', **kw),
	'hpc_local':         lambda X, **kw: run_bnlearn_local(X, method='hpc', **kw),
	'iamb_local':        lambda X, **kw: run_bnlearn_local(X, method='iamb', **kw),
	'iamb_fdr_local':    lambda X, **kw: run_bnlearn_local(X, method='iamb.fdr', **kw),
	'inter_iamb_local':  lambda X, **kw: run_bnlearn_local(X, method='inter.iamb', **kw),
	'mmpc_local':        lambda X, **kw: run_bnlearn_local(X, method='mmpc', **kw),
	'pc_stable_local':   lambda X, **kw: run_bnlearn_local(X, method='pc.stable', **kw),
	'si_hiton_pc_local': lambda X, **kw: run_bnlearn_local(X, method='si.hiton.pc', **kw),
}


def run_method(method_name, X, **kwargs):
	if method_name not in all_methods:
		raise ValueError(f'Unsupported method: {method_name}')
	return all_methods[method_name](X, **kwargs)
