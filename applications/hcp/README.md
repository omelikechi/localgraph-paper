# Human connectome data analysis

This folder contains code and data for reproducing the human connectome and cognition analysis in [*Local graph estimation: Interpretable network discovery for complex data*](https://github.com/omelikechi/localgraph-paper).

---

## Contents

- `data/`: Raw and cleaned data files, scripts for data cleaning, files for enrichment analysis
- `enrichment/`: Scripts for running enrichment analyses
- `figures/`: Scripts for plotting local graphs
- `results/`: Outputs from `run_pfs.py` and `run_methods.py`
- `run_methods.py`: Script for applying other graph estimation methods to the data
- `run_pfs.py`: Script for applying PFS to the data

---

## Notes
- Data were obtained from the Human Connectome Project (HCP Young Adult cohort) via [ConnectomeDB](https://www.humanconnectome.org), accessed through the [BALSA data portal](https://balsa.wustl.edu)
- PFS is a stochastic algorithm, and results may vary slightly across different machines even when using the same random seed. On a single machine, results are deterministic across runs with a fixed seed.


