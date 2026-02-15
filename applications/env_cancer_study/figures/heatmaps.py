# Plot heatmaps for target variables and select environmental exposure and social variables (Figure 3 in the paper)

import os
import geopandas as gpd
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#--------------------------------
# Configuration
#--------------------------------
save_fig = False
dpi = 300

# specify feature types ('targets', 'exposures', or 'social')
feature_type = 'targets'

if feature_type == 'targets':
	features_to_plot = ['Incidence', 'Mortality']
	nrows = 1
elif feature_type == 'exposures':
	features_to_plot = ['pm2.5(A)', 'Hg(W)', 'SO2(A)', 'TCE(A)']
	nrows = 2
elif feature_type == 'social':
	features_to_plot = ['Poverty', 'Education', 'Hispanic', 'Smoking']
	nrows = 2

fig_name = f'heatmaps_{feature_type}'
ncols = 2
figsize = (16,9)

# Color map
cmap_base = plt.colormaps['Spectral'].reversed()
cmap = mcolors.ListedColormap(cmap_base(np.linspace(1/4, 1, 256)))

# format feature names
feature_renames = {
	'Incidence': 'Cancer incidence',
	'Mortality': 'Cancer mortality',
	'pm2.5(A)': 'PM$\\bf{_{2.5}}$',
	'Hg(W)': 'Mercury (Hg)',
	'SO2(A)': 'Sulfur dioxide (SO$\\bf{_{2}}$)',
	'SO4(W)': 'Sulfate (SO$\\bf{_{4}}$)',
	'TCE(A)': 'TCE (C$\\bf{_{2}}$HCl$\\bf{_{3}}$)',
	'PSATest': 'PSA test',
	'PapSmear': 'Pap smear'
}

#--------------------------------
# Plotting function
#--------------------------------
def plot_heatmaps(
	feature_names,
	file,
	base_path="../data/raw_data",
	shapefile_path="./utils/tl_2022_us_county.shp",
	nrows=2,
	ncols=2,
	figsize=(16,9),
	cmap=cmap,
	fig_name=None,
	save_fig=False,
	dpi=300
):
	# Load feature name mapping
	mapping_path = os.path.join(base_path, file + "_feature_names.xlsx")
	mapping_df = pd.read_excel(mapping_path, engine="openpyxl", dtype=str)
	reverse_map = dict(zip(
		mapping_df["Updated Variable Name"].str.strip(),
		mapping_df["Variable Name"].str.strip()
	))

	# Load county shapefile
	gdf = gpd.read_file(shapefile_path)
	gdf["FIPS"] = (gdf["STATEFP"] + gdf["COUNTYFP"]).astype(str).str.zfill(5)

	# Fix outdated Connecticut FIPS codes
	ct_fips_fix = {
		"09110": "09001", "09120": "09003", "09130": "09005", "09140": "09007",
		"09150": "09009", "09160": "09011", "09170": "09013", "09180": "09015"
	}
	gdf["FIPS"] = gdf["FIPS"].replace(ct_fips_fix)
	gdf = gdf[~gdf["STATEFP"].isin(["02", "15"])]  # remove Alaska and Hawaii

	# Load environmental data
	data_path = os.path.join(base_path, file + ".csv")
	df = pd.read_csv(data_path, dtype={'FIPS': str})
	df['FIPS'] = df['FIPS'].str.zfill(5)

	# Setup figure
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
	axes = axes.flatten()

	# Plot each feature
	for i, feature_name in enumerate(feature_names):
		feature_col = reverse_map.get(feature_name, feature_name)
		plot_col = feature_col

		df_feature = df[['FIPS', plot_col]].dropna()
		merged = gdf.merge(df_feature, on="FIPS", how="left").dropna(subset=[plot_col])
		merged[plot_col] = pd.to_numeric(merged[plot_col], errors='coerce')

		# Scale color range
		vmin = merged[plot_col].quantile(0.05)
		vmax = merged[plot_col].quantile(0.95)
		norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

		# Use reversed colormap for Education
		current_cmap = cmap.reversed() if feature_name == 'Education' else cmap

		ax = axes[i]
		merged.plot(column=plot_col, cmap=current_cmap, linewidth=0.3,
					edgecolor="gray", ax=ax, legend=False, norm=norm)

		# Set title
		title = feature_renames.get(feature_name, feature_name)
		ax.set_title(title, fontsize=26, fontweight='bold', pad=1)

		# Standardize layout
		ax.set_xlim(merged.total_bounds[0], merged.total_bounds[2])
		ax.set_ylim(23, 50)
		ax.axis("off")
		ax.set_aspect(1.25)

	# Turn off any unused subplots
	for j in range(i + 1, len(axes)):
		axes[j].axis("off")

	wspace = 0.05 if nrows == 1 else -0.03
	plt.subplots_adjust(left=0, right=1, top=0.95, bottom=0, wspace=wspace, hspace=0.075)

	if save_fig:
		plt.savefig(f'{fig_name}_dpi{dpi}.png', dpi=dpi)
	plt.show()

#--------------------------------
# Call the function
#--------------------------------
plot_heatmaps(
	feature_names=features_to_plot,
	file="eqi2000",
	base_path="../data/raw_data",
	shapefile_path="./utils/tl_2022_us_county.shp",
	ncols=ncols,
	nrows=nrows,
	figsize=figsize,
	fig_name=fig_name,
	save_fig=save_fig,
	dpi=dpi
)
