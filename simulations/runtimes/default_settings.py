# Default settings for runtime experiments

default_settings = {}

# simulation params
default_settings['n'] = 500

# method args
default_settings['pfs_l1_args'] = {'selector':'l1', 'qpath_max':0.25, 'max_radius':3}
default_settings['pfs_gb_args'] = {'selector':'gb', 'qpath_max':0.25, 'max_radius':3}
default_settings['bnlearn_args'] = {}
default_settings['bnlearn_local_args'] = {'radius':3}
default_settings['huge_args'] = {}
default_settings['silggm_args'] = {}