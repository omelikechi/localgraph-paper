# Local graph estimation – paper repository

This repository accompanies the paper:  

**Local graph estimation: Interpretable network discovery for complex data**  
[[arXiv]](https://doi.org/10.48550/arXiv.2507.17172)

It contains all code and data needed to reproduce the analyses, figures, and results presented in the main text and supplementary materials.

> **Note:** This repository is intended to reproduce results in the paper. For a general-use Python package implementing local graph estimation and pathwise feature selection (PFS), see: https://github.com/omelikechi/localgraph

---

## Quickstart

**Simulation illustration** (no data download required, runs in ~30 seconds):
```bash
cd simulations
python illustration.py
```

**HCP neuroimaging application** (~5 minutes):
```bash
cd applications/hcp
python run_pfs.py
```

All scripts must be run from within their own directory as shown above.

---

## Installation

Python 3.10+ is required. Install all dependencies with:
```bash
pip install -r requirements.txt
```

- Some comparison methods require a working R installation.
- If R-based methods are used, ensure R is installed and accessible from your system path.

---

## Repository structure

- `applications/`: peproducible analyses for all real-data applications:
	- `alzheimers/`: Alzheimer's single-nucleus RNA-seq
	- `breast_cancer/`: TCGA breast cancer multi-omics
	- `env_cancer_study/`: environmental cancer study
	- `hcp/`: HCP cognition study

	Each application folder contains:
	- `run_pfs.py`: runs PFS on the dataset
	- `run_methods.py`: runs comparison methods
	- `data/`: preprocessing code and cleaned data
	- `figures/` and `enrichment/`: plotting and enrichment scripts where applicable

- `simulations/`: code for simulation studies and the illustrative example
- `methods/`: wrappers and helper scripts for comparison methods
- `additional_figures/`: script for the conceptual illustration in the supplement
- `requirements.txt`: python dependencies

---

## Datasets

Cleaned datasets are included for all applications except the Alzheimer's disease dataset, which exceeds GitHub's file size limits. Links to all original datasets are provided in the paper.

---

## Contact

For questions about the code or method, please contact:  
Omar Melikechi at omar.melikechi@duke.edu


