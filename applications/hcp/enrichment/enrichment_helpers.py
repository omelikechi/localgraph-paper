# Helper functions for HCP enrichment analyses

from collections import Counter
from scipy.stats import hypergeom

#----------------------------------------------------------------
# DK feature identification and normalization
#----------------------------------------------------------------
def is_cortical_dk_feature(name):
	return (
		(name.startswith('FS_L_') or name.startswith('FS_R_')) and
		(name.endswith('_Thck') or name.endswith('_Area'))
	)


def normalize_fs_name(name):
	"""
	FS_L_Caudalmiddlefrontal_Thck -> 'Left caudalmiddlefrontal'
	"""
	hemi = 'Left' if name.startswith('FS_L_') else 'Right'
	core = name.replace('FS_L_', '').replace('FS_R_', '')
	core = core.replace('_Thck', '').replace('_Area', '')
	return f'{hemi} {core.lower()}'


#----------------------------------------------------------------
# ROI -> Yeo system mapping
#----------------------------------------------------------------
def build_roi_to_system(feature_names, dk_to_yeo):
	"""
	Build mapping from ROI name (e.g. 'Left caudalmiddlefrontal')
	to Yeo-7 system label.
	"""
	roi_to_system = {}

	for fname in feature_names:
		if not is_cortical_dk_feature(fname):
			continue

		roi = normalize_fs_name(fname)
		key = roi.split(' ', 1)[1]
		roi_to_system[roi] = dk_to_yeo.get(key, 'Unknown')

	return roi_to_system


#----------------------------------------------------------------
# Hypergeometric over-representation test for Yeo-7 systems
#----------------------------------------------------------------
def yeo_enrichment(nodes, roi_to_system):
	nodes = list(set(nodes))

	background = Counter(roi_to_system.values())
	foreground = Counter(
		roi_to_system[n] for n in nodes if n in roi_to_system
	)

	N = len(roi_to_system)
	n = len(nodes)

	results = {}
	for system, k in foreground.items():
		K = background[system]
		pval = hypergeom.sf(k - 1, N, K, n)
		results[system] = {'overlap': k, 'total': K, 'p_value': pval}

	return results


