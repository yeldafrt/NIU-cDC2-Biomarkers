# Explainable Machine Learning-Based Identification of Transcriptomic Biomarkers in CD1c+ Dendritic Cells for Non-Infectious Uveitis

This repository contains the complete analysis package for the study: **"Explainable Machine Learning-Based Identification of Transcriptomic Biomarkers in CD1c+ Dendritic Cells for Non-Infectious Uveitis: An Integrative Analysis of Bulk RNA-Seq Data"**.

## Overview

An explainable machine learning framework was applied to bulk RNA-seq data from CD1c+ conventional dendritic cells type 2 (cDC2) isolated from peripheral blood of non-infectious uveitis (NIU) patients and healthy controls. Two independent cohorts (GSE194060 and GSE195501; 78 samples total) were integrated under a strict leakage-free pipeline to identify a parsimonious 20-gene diagnostic signature.

**Key results:**

- Cross-validation AUC: 0.916 | Independent test AUC: 0.855
- Permutation test: p = 0.005
- LODO cross-platform validation: mean AUC = 0.82
- All 20 selected genes are disease-upregulated
- Top biomarkers: CD180, CX3CR1, CCR2, TLR7, DHRS9
- 235 significantly enriched pathway terms (adjusted p < 0.05)

## Repository Structure

```
.
├── code/
│   ├── 01_preprocessing/       # Data download, normalization, batch correction
│   ├── 02_model_training/      # Feature selection, model training, grid search
│   ├── 03_evaluation/          # Test evaluation, bootstrap, permutation test, LODO
│   ├── 04_enrichment/          # Gene annotation and pathway enrichment analysis
│   └── 05_figures/             # Figure generation scripts
├── csv/
│   ├── grid_search_results.csv
│   ├── model_comparison.csv
│   ├── ablation_results.csv
│   ├── lodo_cv_results.csv
│   ├── subgroup_analysis.csv
│   ├── threshold_analysis.csv
│   ├── shap_feature_importance.csv
│   ├── permutation_importance_selected.csv
│   ├── train_metadata.csv
│   ├── test_metadata.csv
│   ├── bootstrap_confidence_intervals.json
│   ├── permutation_test.json
│   ├── overfitting_analysis.json
│   └── final_summary.json
├── dataset/
│   ├── raw_data.zip            # Raw count matrices from GEO
│   └── preprocessed_datasets.zip
│       ├── train_full_8767genes.csv
│       ├── test_full_8767genes.csv
│       ├── train_selected_20genes_scaled.csv
│       ├── test_selected_20genes_scaled.csv
│       ├── train_metadata.csv
│       ├── test_metadata.csv
│       ├── gene_list_8767.csv
│       └── gene_list_selected_20.csv
├── figures/
│   ├── Figure_1_package_v61.zip    # Pipeline overview diagram
│   ├── Figure_2_v61_package.zip    # Permutation test + feature ablation
│   └── Figure_3_v61_package.zip    # SHAP + permutation importance
├── gene_enrichment_results/
│   ├── enrichment/             # Per-library enrichment results
│   ├── enrichment_all_combined.csv
│   ├── enrichment_significant_only.csv
│   ├── gene_annotation_full.csv
│   ├── genes_20_annotated_with_importance.csv
│   ├── table_for_paper_20genes.csv
│   ├── gene_list_all20.txt
│   └── gene_list_disease_up.txt
├── models/
│   ├── BEST_MODEL_v6.joblib    # Final L2-LR model (scikit-learn Pipeline)
│   ├── BEST_MODEL_v6.pkl
│   ├── LR_L2.joblib
│   ├── SVM_Linear.joblib
│   ├── SVM_RBF.joblib
│   └── Random_Forest.joblib
├── npy/
│   ├── X_train.npy             # Training set (62 x 8767)
│   ├── X_test.npy              # Test set (16 x 8767)
│   ├── X_train_scaled.npy      # Scaled training set (62 x 20)
│   ├── X_test_scaled.npy       # Scaled test set (16 x 20)
│   ├── feature_names.npy       # 8767 ENSEMBL gene IDs
│   ├── selected_features.npy   # 20 selected ENSEMBL gene IDs
│   ├── shap_values_train.npy
│   ├── shap_values_test.npy
│   ├── y_pred_test.npy
│   ├── y_pred_train.npy
│   ├── y_proba_test.npy
│   └── y_proba_train.npy
├── plot_data/
│   ├── permutation_scores.npy
│   ├── learning_curve.json
│   └── shap_top30.json
└── logs/
    └── pipeline.log
```

## Data Sources

| Dataset | Samples | Platform | Source |
|:---|:---:|:---|:---|
| GSE194060 | 42 (28 NIU + 14 HC) | Illumina NovaSeq 6000 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194060) |
| GSE195501 | 36 (23 NIU + 13 HC) | Illumina NextSeq 550 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195501) |

Original data from: Hiddingh et al. (2023), *eLife*, 12, e74913.

## Preprocessing Pipeline

1. Raw counts merged across 78 samples (65,217 genes)
2. CPM normalization followed by log2(CPM+1) transformation
3. OLS regression-based batch correction (fitted on training set only)
4. Variance-based gene filtering to top 8,767 genes (fitted on training set only)
5. 80/20 batch-aware stratified split (train: 62, test: 16)
6. ANOVA F-test feature selection to 20 genes (fitted on training set only)
7. StandardScaler normalization (fitted on training set only)

## Selected 20-Gene Signature

| Gene Symbol | ENSEMBL ID | SHAP Importance | Direction |
|:---|:---|:---:|:---|
| CD180 | ENSG00000134061 | 0.0148 | Disease-associated |
| CX3CR1 | ENSG00000168329 | 0.0133 | Disease-associated |
| CCR2 | ENSG00000121807 | 0.0128 | Disease-associated |
| TLR7 | ENSG00000196664 | 0.0127 | Disease-associated |
| DHRS9 | ENSG00000073737 | 0.0113 | Disease-associated |
| TLR10 | ENSG00000174123 | 0.0113 | Disease-associated |
| CCR5 | ENSG00000160791 | 0.0112 | Disease-associated |
| ADRB2 | ENSG00000169252 | 0.0111 | Disease-associated |
| MED31 | ENSG00000108590 | 0.0110 | Disease-associated |
| QNG1 | ENSG00000165118 | 0.0109 | Disease-associated |
| RGS18 | ENSG00000150681 | 0.0107 | Disease-associated |
| GAPT | ENSG00000175857 | 0.0101 | Disease-associated |
| ZNF761 | ENSG00000160336 | 0.0095 | Disease-associated |
| KRCC1 | ENSG00000172086 | 0.0093 | Disease-associated |
| GIMAP2 | ENSG00000106560 | 0.0092 | Disease-associated |
| VTA1 | ENSG00000009844 | 0.0089 | Disease-associated |
| ERGIC2 | ENSG00000087502 | 0.0086 | Disease-associated |
| CD200R1 | ENSG00000163606 | 0.0085 | Disease-associated |
| TLR8 | ENSG00000101916 | 0.0084 | Disease-associated |
| CYSLTR1 | ENSG00000173198 | 0.0081 | Disease-associated |

## Requirements

- Python 3.11+
- scikit-learn
- numpy, pandas
- matplotlib, seaborn
- shap
- gseapy
- mygene
- joblib

## Usage

To load the trained model and make predictions:

```python
import joblib
import numpy as np

model = joblib.load("models/BEST_MODEL_v6.joblib")
X_test = np.load("npy/X_test.npy")

predictions = model.predict(X_test)
probabilities = model.predict_proba(X_test)[:, 1]
```

## Citation

If you use this code or data, please cite the associated publication (details to be added upon acceptance).
