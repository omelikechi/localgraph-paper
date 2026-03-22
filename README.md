# Local graph estimation – paper repository

This repository accompanies the paper:  
**Local graph estimation: Interpretable network discovery for complex data**  
[[arXiv]](https://doi.org/10.48550/arXiv.2507.17172)

It contains all code and data needed to reproduce the analyses, figures, and results presented in the main text and supplementary materials.

> **Note:** This repository is intended to reproduce results in the paper. For a general-use Python package implementing local graph estimation and pathwise feature selection (PFS), see: https://github.com/omelikechi/localgraph

---

## Quickstart

To verify the code works, we recommend starting with the simulation illustration and one real-data application.

**Simulation illustration** (no data download required, runs in ~20 seconds):
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
Omar Melikechi — omar.melikechi@duke.edu




















<!-- # Local graph estimation – paper repository

This repository accompanies the paper:

**Local graph estimation: Interpretable network discovery for complex data**<br>
[[arXiv]](https://doi.org/10.48550/arXiv.2507.17172)

It contains all code and data needed to reproduce the analyses, figures, and results presented in the main text and supplementary materials.

> **Note:** This repository is intended to reproduce results in the paper. For a general-use Python package implementing local graph estimation and pathwise feature selection (PFS), see:  
https://github.com/omelikechi/localgraph

---

## Repository structure

- applications/  
	Reproducible analyses for all real-data applications in the paper:
	- `alzheimers/`
	- `breast_cancer/`
	- `env_cancer_study/`
	- `hcp/`

	Each application folder contains:
	- `run_pfs.py` and `run_methods.py` scripts
	- preprocessing code in `data/`
	- enrichment and plotting scripts where applicable

- simulations/  
	Code for both simulation studies and the illustrative example used in the paper.

- methods/  
	Wrappers and helper scripts used to run comparison methods.

- additional_figures/  
	Script for generating the conceptual illustration in the supplement.

- requirements.txt  
	Python dependencies required to reproduce all experiments.

---

## Installation

Install all required dependencies with:

```
pip install -r requirements.txt
```

Note:
- Some comparison methods require a working R installation.
- If R-based methods are used, ensure R is installed and accessible from your system path.

---

## Datasets

Cleaned datasets used for analysis are available for each application except the Alzheimer's dataset, which exceeds GitHub's capacity. Links to the original datasets associated with each application are provided in the paper.

---

## Questions

For issues related to the method or repository, please contact:

Omar Melikechi  
omar.melikechi@duke.edu


 -->