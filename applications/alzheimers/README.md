# Alzheimer's disease data analysis

This folder contains code and data for reproducing the Alzheimer's disease (AD) analyses (Section 7).

---

## Contents

- `data/`: Script for data cleaning (data too large for GitHub; see not below)
- `enrichment/`: Scripts for running gene enrichment analyses
- `figures/`: Scripts for plotting local graphs
- `results/`: Outputs from `run_pfs.py` and `run_methods.py`
- `run_methods.py`: Script for applying other graph estimation methods to the data
- `run_pfs.py`: Script for applying PFS to the data

---

## Notes
- The raw and cleaned data files exceed GitHub’s 100MB file limit and therefore cannot be included in the public repo. However, these data can be downloaded for free from [https://adsn.ddnetbio.com](https://adsn.ddnetbio.com).
- PFS is a stochastic algorithm, and results may vary slightly across different machines even when using the same random seed. On a single machine, results are deterministic across runs with a fixed seed.


