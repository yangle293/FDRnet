The code in this directory is preprocessing scripts for our FDRnet paper.

1. Breast Cancer study: A pre-processing pipeline for computing significance scores for the copy number data and combining scores derived from mutation and copy number data. This pipeline is used to generate input data (local FDR scores) for FDRnet in breast cancer study.

2. Lymphoma study: A R script for computing local FDR scores from p-values that are obtained from gene differential expression analysis. We use R package "locfdr" inside.
