# Environmental and social drivers of cancer

This folder contains code and data for reproducing the county-level cancer analyses (Section 4).

---

## Contents

- `data/`: Datasets and scripts for data cleaning
- `figures/`: Scripts for plotting local graphs
- `results/`: Outputs from `run_pfs.py` and `run_methods.py`
- `run_methods.py`: Script for applying other graph estimation methods to the data
- `run_pfs.py`: Script for applying PFS to the data

---

## Notes
- **Sources:**
	- **Cancer outcomes and screening:** [State Cancer Profiles](https://statecancerprofiles.cancer.gov/)
	- **Demographics:** U.S. Census Bureau (2000 decennial census)
	- **Environmental exposures:** [EPA Environmental Quality Index (EQI)](https://www.epa.gov/enviroatlas/environmental-quality-index)
- PFS is a stochastic algorithm, and results may vary slightly across different machines even when using the same random seed. On a single machine, results are deterministic across runs with a fixed seed.


