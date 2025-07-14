# Environmental and social drivers of cancer

This folder contains code and data for reproducing the county-level cancer analysis in Section 3 of [*Local graph estimation: Interpretable network discovery for complex data*](https://github.com/omelikechi/localgraph-paper). The analysis applies pathwise feature selection (PFS) to a combination of environmental exposures, socioeconomic factors, and demographic characteristics to study county-level cancer incidence and mortality rates across the contiguous United States.

---

## Contents

- `run_pfs.py`: Main script for applying PFS to the county-level data  
- `figures/`: Scripts for generating all figures, including:
	- `plot_pfs_results.py`: Reproduces **Figure 2** (PFS estimate of the local graph)
	- `plot_glasso_results.py`: Reproduces **Figure 1e** (Graphical lasso applied to the data)
	- `heatmaps.py`: Reproduces the heatmaps in **Figure 3**
	- `scatter_plots.py`: Reproduces **Figure S2** in the Supplement
	- `utils/`: Includes shapefiles for plotting U.S. county heatmaps
- `data/`: Contains raw and cleaned data
- `results/`: Contains output from `run_pfs.py`, including the final estimated graph

---

## Data Overview

**Sources:**
- **Cancer outcomes and screening:** [State Cancer Profiles](https://statecancerprofiles.cancer.gov/)
- **Demographics:** U.S. Census Bureau (2000 decennial census)
- **Environmental exposures:** [EPA Environmental Quality Index (EQI)](https://www.epa.gov/enviroatlas/environmental-quality-index)

### Raw data:
- Cancer incidence and mortality (2017–2022)
- Cancer screening and smoking rates (2017–2019)
- County-level environmental and social indicators (2000–2005)
- County-level race/ethnicity, sex, and age distributions

### Cleaned data (`data/cleaned_data/`)
- `cleaned_data.csv`: Final dataset including the two target variables (cancer incidence and mortality) followed by all filtered covariates. Column names are preserved in the header row.

**Filtering and cleaning:**
- Excluded Indiana and Kansas because they did not report incidence
- Excluded Union County, FL due to extreme incidence and mortality values
- Dropped variables missing from ≥ 5 states or ≥ 1500 counties
- Final dataset: 2857 counties × 165 variables

---

## PFS Implementation

- Base selector: Random forests (scikit-learn)  
- Maximum radius: `2`  
- Pathwise FDR threshold: `1.0`  
- Neighborhood FDR thresholds:
	- Incidence: `0.04`
	- Mortality: `0.009`
	- Hispanic: `0.01`
	- Poverty: `0.011`
	- Smoking: `0.01`
	- Education: custom threshold using efp scores for tie-breaking
	- All others: `0.008`

---

## Results

- `env_study_pfs_results.pkl`: Main output from `run_pfs.py`
- `adjusted_graph_pfs.graphml`: PFS graph with node positions adjusted for better visibility (used for **Figure 3**)
- `env_study_glasso_results.pkl`: Results from applying the graphical lasso to the data
- `adjusted_graph_glasso.graphml`: Graphical lasso graph with node positions adjusted (used for **Figure 1e**)

---

## Notes

PFS is a stochastic algorithm, and results may vary slightly across different machines even when using the same random seed. On a single machine, results are deterministic across runs with a fixed seed.
