# Format the HCP feature names for cleaner plotting and analysis

def format_feature_names(var):
	# drop age-adjustment suffix
	var = var.replace('_AgeAdj', '')

	# target variable
	if var == 'CogFluidComp':
		return 'Fluid\nCognition'

	# behavioral / phenotype summaries
	simple_map = {
		'WM_Task_Acc': 'Working\nmemory',
		'Language_Task_Acc': 'Language\ntask',
		'Emotion_Task_Acc': 'Emotion\ntask',
		'Endurance': 'Endurance',
		'Dexterity': 'Dexterity',
		'Strength': 'Strength',
		'GaitSpeed_Comp': 'Gait Speed',
		'Noise_Comp': 'Noise\nSensitivity',
		'Odor': 'Olfaction',
		'Taste': 'Taste',
		'Mars_Log_Score': 'Visual\nAcuity',
		'NEOFAC_A': 'Agreeableness',
		'NEOFAC_C': 'Conscientiousness',
		'NEOFAC_E': 'Extraversion',
		'NEOFAC_N': 'Neuroticism',
		'NEOFAC_O': 'Openness',
	}
	if var in simple_map:
		return simple_map[var]

	# FreeSurfer variables
	if var.startswith('FS_'):
		v = var[3:]  # drop FS_

		hemi = None
		if v.startswith('L_'):
			hemi = 'Left\n'
			v = v[2:]
		elif v.startswith('R_'):
			hemi = 'Right\n'
			v = v[2:]

		parts = v.split('_')

		measure = parts[-1]
		region = ' '.join(parts[:-1])

		measure_map = {
			'Thck': '(thick)',
			'Area': '(area)',
			'Vol': '(vol)'
		}
		measure = measure_map.get(measure, measure)

		if hemi:
			return f'{hemi}{region}\n{measure}'
		else:
			return f'{region}\n{measure}'

	# fallback
	return var