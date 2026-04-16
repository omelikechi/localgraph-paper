# TCGA breast cancer analysis

This folder contains code and data for reproducing the breast cancer analyses (Section 5).

---

## Contents

- `data/`: Datasets and scripts for data cleaning
- `enrichment/`: Scripts for running gene and protein enrichment analyses
- `figures/`: Scripts for plotting local graphs
- `results/`: Outputs from `run_pfs.py` and `run_methods.py`
- `run_methods.py`: Script for applying other graph estimation methods to the data
- `run_pfs.py`: Script for applying PFS to the data

---

## Notes
- The RNA-seq data file exceeds GitHub’s 100MB file limit and is not included in the public repo. However, it can be downloaded for free from [LinkedOmics](https://www.linkedomics.org/data_download/TCGA-BRCA/).
- PFS is a stochastic algorithm, and results may vary slightly across different machines even when using the same random seed. On a single machine, results are deterministic across runs with a fixed seed.

