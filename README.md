# Anthrax-MLNet: A Hybrid Machine Learning Approach with Multi-Step Feature Selection for Identifying Candidate Biomarkers and Drug Targets in Anthrax

This repository contains the code, processed datasets, and figures for the paper **"Anthrax-MLNet: A Hybrid Machine Learning Approach with Multi-Step Feature Selection for Identifying Candidate Biomarkers and Drug Targets in Anthrax"**.

> **Note:** This is a small-sample, computational discovery study (n = 11). All findings are **candidate** biomarkers and drug targets. They require validation in larger independent cohorts and experimental follow-up before any diagnostic or therapeutic claim can be made.

---

## Overview

Anthrax, caused by *Bacillus anthracis*, progresses rapidly and remains a concern for infectious disease management and biodefense. Transcriptomic datasets for this pathogen are small, which limits both biomarker discovery and validation.

This project integrates two public GEO series and applies a leakage-free, multi-step analysis pipeline:

| Dataset | Description | Samples |
| --- | --- | --- |
| **GSE14390** | Spore-exposed alveolar macrophages | 7 |
| **GSE12131** | Lethal-toxin-treated endothelial cells | 4 |
| **Total** | 5 anthrax-exposed, 6 control | **11** |

## Pipeline

1. **Preprocessing and harmonization:** shared genes retained, probes averaged per gene symbol, GSE12131 log2-transformed, per-cohort gene-wise z-scaling before pooling.
2. **Differential expression:** Welch's t-test with Benjamini–Hochberg FDR (adj. P < 0.05, |logFC| > 0.5), complemented by a limma workflow.
3. **Multi-step feature selection:** LASSO logistic regression, Random Forest importance, and Recursive Feature Elimination (RFE), combined with the DEG list by consensus. Selection is repeated inside each training fold to avoid leakage. The result is a **12-gene candidate panel**.
4. **Functional enrichment:** KEGG, WikiPathways, Reactome, BioPlanet, and Gene Ontology.
5. **PPI network and hub genes:** STRING via NetworkAnalyst, Cytoscape (v3.10.3), and CytoHubba (degree centrality).
6. **Regulatory networks:** TF–gene (ENCODE), TF–miRNA co-regulation (RegNetwork), and gene–miRNA (miRTarBase v9.0).
7. **Internal validation:** leave-one-out cross-validation (L2 logistic regression, standardization and fitting confined to training folds), 5000-label permutation test, repeated 3-fold CV, and out-of-bag bootstrap.
8. **Effect size and power:** Cohen's d, achieved power, and required sample size per gene.
9. **Bootstrap clustering stability:** multiscale bootstrap (pvclust, 1000 replicates) and MDS over 500 gene-resampling replicates.
10. **Drug signatures and docking:** DSigDB via Enrichr, then AutoDock Vina docking of 5 compounds against 5 hub proteins.

## Key Results

### Classification performance (LOOCV, n = 11)

| Metric | 12 overlapping DEGs | 5 hub genes |
| --- | --- | --- |
| AUC (DeLong 95% CI) | **0.967** (0.874–1.000) | 0.667 (0.305–1.000) |
| Accuracy | 90.9% | 72.7% |
| Sensitivity / Specificity | 100% / 83.3% | 100% / 50% |
| Permutation P | **0.003** | 0.216 |
| Majority-class baseline | 54.5% | 54.5% |

- The 12-gene AUC was reproduced by repeated 3-fold CV (0.968) and out-of-bag bootstrap (0.967).
- **GULP1, GATAD2A, and ACSS1** showed the highest selection stability across folds (11, 11, and 9 folds respectively).
- Bootstrap clustering split samples by **infection status** (agreement 0.82) rather than **cohort of origin** (0.55), so the signal is not simply a batch artifact.
- 10 of the 12 panel genes shifted in the same direction in both cohorts.

### Hub genes vs. diagnostic value

- Top hub genes by degree: **WNT6, FZD9, DKK4, FZD6, FZD3**, all converging on **Wnt signaling**.
- As a group they showed **no discriminative capacity** (AUC 0.667, P = 0.22). Topological centrality and diagnostic value are distinct properties.
- WNT6 is the exception because it is also in the DEG list (d = 2.07, power 0.86).
- The hub genes are kept as structurally central nodes and as docking receptors, not as candidate markers.

### Well-powered core of the panel

LOC101928461, LINC00485, WNT6, RNF128, and LOC100506530 reached FDR < 0.05. VWC2 reached adequate power but not FDR < 0.05. Other genes (e.g., ACSS1, IL20RA, GATAD2A, TTC6) were underpowered at this sample size.

### Drug signatures and molecular docking

- Candidate compounds from DSigDB: **irinotecan, daunorubicin, celastrol, hesperidin, quercetin**.
- Predicted binding energies ranged from **−6.8 to −11.1 kcal/mol**. 24 of 25 pairs were below −7.0 kcal/mol.
- Strongest pair: **irinotecan–FZD6 (−11.1 kcal/mol)**. Irinotecan had the best mean affinity (−9.72 kcal/mol).
- Docking used a grid centered on each receptor, not a validated binding site, so these are **non-site-specific, preliminary** results that need MD simulation and experimental assays.

## Repository Structure

Anthrax_Analysis/
├── Analysis_Code/       # Scripts for preprocessing, DEG analysis, feature selection,
│                        # validation, bootstrap, effect size/power analysis
├── Analysis_Datasets/   # Processed expression matrices, DEG lists, selected biomarker panel
├── Analysis_Figure/     # Figures generated for the manuscript (volcano, heatmap,
│                        # stability, enrichment, PPI, ROC, bootstrap, docking)
├── Supplementary File/  # Supplementary tables and figures referenced in the paper
│                        # (e.g., Tables S1–S12, Figures S1–S6)
└── README.md

## Requirements

The analysis was carried out with a combination of Python and R, plus web-based and desktop tools.

**Python:** `numpy`, `pandas`, `scikit-learn`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`

**R:** `limma`, `pvclust`, `pROC`

**External tools and databases:**
NCBI GEO, NetworkAnalyst, STRING, Cytoscape (v3.10.3) with CytoHubba, Enrichr (DSigDB), ENCODE, RegNetwork, miRTarBase v9.0, RCSB PDB, PubChem, AutoDockTools 1.5.7, AutoDock Vina, PyMOL.

## How to Reproduce

1. Clone the repository:
   ```bash
   git clone https://github.com/MostaryKhatun/Anthrax_Analysis.git
   cd Anthrax_Analysis
   ```
2. Install the required Python and R packages listed above.
3. Download the raw data from GEO (**GSE14390** and **GSE12131**), or use the processed matrices in `Analysis_Datasets/`.
4. Run the scripts in `Analysis_Code/` in order: preprocessing, DEG analysis, feature selection, validation, bootstrap clustering, and effect size/power analysis.
5. Network, enrichment, and drug-signature steps were done on the web platforms listed above. Docking was performed in AutoDock Vina.

## Limitations

- Only 11 samples from two cell-culture models (alveolar macrophages and endothelial cells), not patient material.
- Validation is internal only. No adequately sized external anthrax cohort was available.
- Wide confidence intervals; many panel genes are underpowered.
- Drug and docking results are exploratory and unvalidated.


