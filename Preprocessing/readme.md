# Preprocessing Scripts
The code in this directory are preprocessing scripts for our FDRnet paper.

- `BreastCancerPreprocessing.m` (Breast Cancer study): Implemented the pre-processing pipeline proposed in the paper for computing significance scores for the copy number data and combining scores derived from mutation and copy number data. This pipeline is used to generate input data (local FDR scores) for FDRnet in breast cancer study.

- `LymphomaPreprocessing.R` (Lymphoma study): A R script for computing local FDR scores from p-values that are obtained from gene differential expression analysis. This code require the R package `locfdr`. We use R script here since the original data is available as an R dataset.
