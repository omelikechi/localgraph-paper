# Plot relationships between exposures using cleaned data

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#--------------------------------
# Load cleaned data
#--------------------------------

save_fig = False
dpi = 300
fig_name = f'eqi_nonlinear_dpi{dpi}'

# Load cleaned full data
data_path = "../data/cleaned_data/cleaned_data.csv"
df = pd.read_csv(data_path)
df = df.dropna()

feature_names = df.columns.tolist()
X = df.to_numpy()

#--------------------------------
# Scatter plots
#--------------------------------

y_vars = ['pm2.5(A)', 'SO2(A)', 'Hg(W)']
shared_vars = ['O3(A)', 'SO4(W)', 'Ca(W)', 'TCE(A)', 'Quinoline(A)', 'Acrylonitrile(A)']

var_names = {
    'Hg(W)': 'Hg',
    'pm2.5(A)': 'PM$_{2.5}$',
    'SO2(A)': 'SO$_{2}$',
    'O3(A)': 'O$_3$',
    'SO4(W)': 'SO$_4$',
    'Ca(W)': 'Ca',
    'TCE(A)': 'TCE',
    'Quinoline(A)': 'Quinoline',
    'Acrylonitrile(A)': 'Acrylonitrile'
}

nrows = len(y_vars)
ncols = int(np.ceil(len(shared_vars)))
fig, ax = plt.subplots(nrows, ncols, figsize=(18,8.5), sharex=False)
ax = ax.flatten()

colors = ['deepskyblue', 'orange', 'limegreen']

for row in range(nrows):
    color = colors[row]
    y_var = y_vars[row]
    y_vals = df[y_var].to_numpy()

    for i, x_var in enumerate(shared_vars):
        x_vals = df[x_var].to_numpy()
        plot_idx = row * ncols + i
        ax[plot_idx].scatter(x_vals, y_vals, s=10, color=color, alpha=0.5)
        ax[plot_idx].set_xticks([])
        ax[plot_idx].set_yticks([])

        if i == 0:
            ax[plot_idx].set_ylabel(var_names[y_var], fontsize=24)
        if row == nrows - 1:
            ax[plot_idx].set_xlabel(var_names[x_var], fontsize=24)

plt.tight_layout()
if save_fig:
    plt.savefig(fig_name + ".png", dpi=dpi)
plt.show()
