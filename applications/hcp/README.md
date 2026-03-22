# Human connectome data analysis

This folder contains code and data for reproducing the human connectome and cognition analyses (Section 6).

---

## Contents

- `data/`: Datasets and scripts for data cleaning
- `enrichment/`: Scripts for running enrichment analyses
- `figures/`: Scripts for plotting local graphs
- `results/`: Outputs from `run_pfs.py` and `run_methods.py`
- `run_methods.py`: Script for applying other graph estimation methods to the data
- `run_pfs.py`: Script for applying PFS to the data

---

## Notes
- Data were obtained from the Human Connectome Project (HCP Young Adult cohort) via [ConnectomeDB](https://www.humanconnectome.org), accessed through the [BALSA data portal](https://balsa.wustl.edu)
- PFS is a stochastic algorithm, and results may vary slightly across different machines even when using the same random seed. On a single machine, results are deterministic across runs with a fixed seed.


