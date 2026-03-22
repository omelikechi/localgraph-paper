# TCGA breast cancer analysis

This folder contains code and data for reproducing the breast cancer analysis in [*Local graph estimation: Interpretable network discovery for complex data*](https://github.com/omelikechi/localgraph-paper). The analysis applies pathwise feature selection (PFS) to a combintion of gene (RNAseq), miRNA, and protein (RPPA) data from The Cancer Genome Atlas (TCGA) to uncover local structures associated with clinical outcomes in breast cancer patients.

---

## Contents

- `run_pfs.py`: Main script for applying PFS to the breast cancer data
- `figures/`: Reproduce figures from the paper (local graph around target features)
- `data/`: Contains raw and cleaned data
- `results/`: Contains outputs from `run_pfs.py` and from other methods

---

## Data Overview

**Source:** [LinkedOmics: TCGA-BRCA](https://www.linkedomics.org/data_download/TCGA-BRCA/)

### Raw data (`data/raw_data/`)
- `rnaseq.txt`: RNAseq gene expression (HiSeq, Log2(RPKM+1))  
- `mirna.txt`: miRNA expression (HiSeq, RPM, Log2(Val+1))  
- `rppa.txt`: Protein expression (RPPA, normalized signal)  
- `clinical.txt`: Clinical metadata (including status, stage, histological type)  

**Note:** `rnaseq.txt` exceeds GitHub’s 100MB file limit and is not included in the public repo but can be downloaded for free from the [LinkedOmics](https://www.linkedomics.org/#/).

### Cleaned data (`data/cleaned_data/`)
- `cleaned_data.csv`: Final dataset including the three target clinical variables (histological type, pathologic stage, status) followed by all filtered features across RNA, miRNA, and protein modalities. Feature names are preserved in the header row.

**Filtering steps:**
- Retained genes (RNAseq) with mean expression or variance above the 75th percentile  
- Retained only samples with complete data across all modalities and clinical targets  
- Removed features with absolute pairwise correlation > 0.999 (kept one from each pair)  
- Removed features or samples with missing values  

**Final dimensions:**
- RNAseq: 9785 genes  
- miRNA: 819 features  
- RPPA: 137 proteins  
- Final feature matrix: 547 samples × 10,741 features  
- Clinical targets: 547 samples × 3 targets

---

<!-- ## PFS Implementation

- Base selector: Gradient boosting (XGBoost)  
- Number of subsamples: `B = 200`  
- Preselection expansion factor: `2.5`  
- Maximum radius: `3`  
- Pathwise FDR threshold: `1.0`  
- Neighborhood FDR thresholds:
  - Histological type: `0.25`
  - Pathologic stage: `0.2`
  - Status: `0.2`
- Local FDR thresholds:
  - Radius 1: `0.023`
  - Radius 2: `0.0175`
- Neighborhood threshold for cross-modal edges: `0.035`

--- -->

## Results

- `breast_cancer_pfs_results.pkl`: Main output from `run_pfs.py`  
- `adjusted_graph.graphml`: Graph with node positions manually adjusted for visibility (used for Figure 4)  

---

## Notes

PFS is a stochastic algorithm, and results may vary slightly across different machines even when using the same random seed. On a single machine, results are deterministic across runs with a fixed seed.
