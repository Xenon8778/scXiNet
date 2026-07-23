# scXiNet: Constructing Gene Regulatory Networks using Chatterjee's Rank Correlation with Single-cell Transcriptomic Data

[![R](https://img.shields.io/badge/R-100%25-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Code accompanying the paper *"Constructing Gene Regulatory Network using Chatterjee's Rank Correlation with Single-cell Transcriptomic Data"* (Gupta, Chaudhuri, Raghuraman, Ni & Cai).

## Abstract

Discovering gene regulatory networks (GRNs) from single-cell RNA sequencing (scRNA-seq) data is critical for understanding cellular function, but existing methods are limited by strong theoretical assumptions or high computational complexity. We introduce a multiple testing framework for inference of GRNs using Chatterjee's rank correlation coefficient, a nonparametric measure of dependence. Our approach overcomes the limitations of traditional methods while offering a transparent, scalable, and computationally efficient alternative to recent black-box machine learning models. We address the challenge of non-independent observations in scRNA-seq with a data-driven algorithm for estimating robust testing cutoffs, and exploit the asymmetric nature of Chatterjee's correlation to propose a new test for directed regulation. Our method consistently outperforms state-of-the-art approaches in recovering true regulatory links from both simulated and real datasets.

## Citation

If you use this code, please cite:

```bibtex
@article{gupta2025constructing,
  title={Constructing Gene Regulatory Network using Chatterjee’s Rank Correlation with Single-cell Transcriptomic Data},
  author={Gupta, Shreyan and Chaudhuri, Anamitra and Raghuraman, Vishnuvasan and Ni, Yang and Cai, James J},
  journal={bioRxiv},
  pages={2025--09},
  year={2025},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Pipeline Overview

The codebase follows this general workflow:

1. **`PreprocessingCodes/`** — QC, normalization, and clustering of raw scRNA-seq data (Seurat v5: filtering, `LogNormalize`, PCA/UMAP, Louvain clustering) to produce the gene expression matrices used downstream.
2. **`SimulationCodes/`** — Generates synthetic scRNA-seq datasets with predefined linear/non-linear gene-gene interactions (Poisson-approximation simulator) and SERGIO-based trajectory data, used for cutoff estimation and benchmarking.
3. **`utils/`** — Core library: Chatterjee's correlation computation (parallelized across cores for large gene panels) and the data-driven FPR-controlled cutoff selection (Algorithms 1 & 2 in the paper).
4. **`AnalysisCodes/`** — Benchmarking against Pearson, Spearman, distance correlatiom, mutual information, GRNBoost2, and pcNet; cutoff transferability studies; real-data GRN construction and validation against STRING.

Typical use: preprocess data → estimate cutoffs on a reference/simulated dataset (`utils/` + `SimulationCodes/`) → apply cutoffs to construct and validate GRNs (`AnalysisCodes/`).

## Folder Contents

| Folder | Description |
|---|---|
| `PreprocessingCodes/` | Scripts for pre-processing scRNA-seq datasets (QC, normalization, clustering). |
| `SimulationCodes/` | Scripts to simulate scRNA-seq datasets with known gene-gene interactions (linear, parabolic, exponential, sinusoidal), used for cutoff estimation and method benchmarking. |
| `AnalysisCodes/` | Scripts used for benchmarking dependency measures and generating the results/figures reported in the paper. |
| `utils/` | Core functions for FPR-based cutoff selection and parallelized computation of Chatterjee's correlation for large gene-expression matrices. |

## Requirements

- R (version used: 4.5.2)
- Key packages:
  - `Seurat` (v5) — preprocessing
  - `doParallel` / `foreach` — parallelized Chatterjee's correlation
  - `XICOR` or custom implementation — Chatterjee's ξ coefficient
  - Standard tidyverse (`dplyr`, `ggplot2`) for analysis/plotting

```r
install.packages(c("Seurat", "doParallel", "foreach", "dplyr", "ggplot2"))
```

## Data Availability

- Real scRNA-seq datasets (PBMC and small intestine epithelial cells) — sourced from GEO; accession numbers listed in Supplementary Table 1 of the paper.
- Simulated datasets (SERGIO) — generator available at [github.com/PayamDiba/SERGIO](https://github.com/PayamDiba/SERGIO).
- Ground-truth protein-protein interaction network — [STRING database](https://string-db.org/) v12.0.

No new data was generated for this study.

## License

This project is licensed under the [MIT License](LICENSE).

## Authors & Contact

Shreyan Gupta, Anamitra Chaudhuri, Vishnuvasan Raghuraman, Yang Ni, James J. Cai (Texas A&M University; UT Austin)

Corresponding author: xenon8778@tamu.edu
