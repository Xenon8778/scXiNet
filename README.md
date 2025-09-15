# Constructing Gene Regulatory Network using Chatterjee’s Rank Correlation with Single-cell Transcriptomic Data

## Abstract
Discovering gene regulatory networks (GRNs) from single-cell RNA sequencing (scRNA-seq) data is critical for understanding cellular function, but existing methods are limited by strong theoretical assumptions or high computational complexity. We introduce a multiple testing framework for inference of GRN by using Chatterjee's rank correlation coefficient, a nonparametric measure of dependence. Our approach overcomes the limitations of traditional methods while offering a transparent, scalable, and computationally efficient alternative to recent black-box machine learning models. Crucially, we address the challenge of non-independent observations in scRNA-seq by developing a data-driven algorithm for estimating robust testing cutoffs. Furthermore, we exploit the asymmetric nature of Chatterjee's correlation to propose a new test for directed regulation, enabling the construction of biologically meaningful and directionally informed GRNs. We demonstrate that our method consistently outperforms state-of-the-art approaches in recovering true regulatory links from both simulated and real datasets, providing a powerful tool for dissecting complex GRNs.

## Folder Contents
- AnalysisCodes - contains scripts used for benchmarking and generating results.
- PreprocessingCodes - contains scripts for pre-processing scRNA-seq datasets.
- SimulationCodes - contains scripts used to simulate scRNA-seq datasets with gene-gene interactions.
- utils - contains functions implemented for FPR based cut-off selection and computing chatterjee's correlation for large matrices using mutliple parallel cores.