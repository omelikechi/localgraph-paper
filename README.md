# Local graph estimation – paper repository

This repository accompanies the paper:

**Local graph estimation: Interpretable network discovery for complex data**

It contains all code, data, and figures needed to reproduce the analyses, figures, and results presented in the main text and supplementary materials.

> ⚠️ **Note:** This repository is intended solely for reproducing the paper.  
> For a general-use Python package implementing local graph estimation and pathwise feature selection (PFS), see:  
> https://github.com/omelikechi/localgraph

---

## Repository structure

- `breast_cancer/`  
  Code and data for analyzing TCGA breast cancer multiomic data.

- `env_cancer_study/`  
  Code and data for the exposomic study linking environmental and social factors to U.S. county-level cancer outcomes.

- `simulations/`  
  Scripts and results for both simulation studies, as well as illustrative examples demonstrating the limitations of global error control.  
  See `simulations/README.md` for full details.

- `additional_figures/`  
  Self-contained script for generating the conceptual illustration in the supplement (Figure S1).

- `requirements.txt`  
  Full list of package dependencies, including frozen versions of `localgraph` and `ipss`.

---

## Getting started

We recommend using a virtual environment to avoid version conflicts and ensure reproducibility.

To install all dependencies, run:

```bash 
pip install -r requirements.txt 
```

---

## Datasets

This repository includes all datasets needed to reproduce the results in the paper. Below is a summary of each dataset, where it is stored in the repository, and a link to the original source (when applicable).

### Breast Cancer Multiomic Data

- **Location in repo**: `breast_cancer/data/`
- **Description**: RNA-seq, miRNA, and RPPA (protein expression) data from breast cancer samples.
- **Original source**: [LinkedOmics: TCGA-BRCA](https://www.linkedomics.org/data_download/TCGA-BRCA/)
- **RNA-seq**: The raw RNA-seq text file is too large to be uploaded to GitHub but can be downloaded from the link above

### County-Level Cancer and Environmental Data

- **Location in repo**: `env_cancer_study/data/`
- **Description**: County-level cancer incidence and mortality rates, environmental exposures, socioeconomic factors, demographic characteristics, and screening data.
- **Original sources**:
  - [State Cancer Profiles (NCI & CDC)](https://statecancerprofiles.cancer.gov/)
  - [EPA Environmental Quality Index (EQI)](https://www.epa.gov/enviroatlas/environmental-quality-index)
  - [U.S. Census Bureau (2000)](https://data.census.gov/cedsci/)

---

## Questions

For issues or questions related to the method or this repository, feel free to open an issue on GitHub or email the corresponding author:

**Omar Melikechi**  
omar.melikechi@gmail.com



