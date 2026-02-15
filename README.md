# Local graph estimation – paper repository

This repository accompanies the paper:

**Local graph estimation: Interpretable network discovery for complex data**

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


