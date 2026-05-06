# Longitudinal Proteomics: Moderated F-test Two-Stage Framework

This repository contains **R code** to reproduce the simulations and analyses for the manuscript:

**Robust Biomarker Discovery and Prediction in Longitudinal Proteomics: A Moderated F-Test–Based Two-Stage Framework**  

## Overview

High-dimensional longitudinal proteomics studies often have **very small numbers of subjects** and repeated measurements over time. Standard approaches (e.g., running a t-test at each time point and picking the smallest p-value) can inflate false positives and fail to capture the joint temporal pattern.

This project implements a **two-stage framework**:

1. **Biomarker discovery (feature selection):**  
   Use **limma**-style linear models with **empirical Bayes moderation** and a **moderated F-test** to test whether a protein shows any treatment-associated change across the entire time course.

2. **Prediction model:**  
   Use selected proteins as features in a **Naive Bayes classifier** to classify treatment status with a **small, practical protein panel**.

## Contents

- **Simulation 1 (cross-sectional):** demonstrates why empirical Bayes moderation is crucial for small sample sizes (moderated t-test vs ordinary t-test).
- **Simulation 2 (longitudinal):** compares ordinary t-test / moderated t-test (min p-value across time) vs moderated F-test (global time-course test), including downstream classification performance.
- **Application:** equine anti-doping longitudinal proteomics study (rhEPO / Mircera).

## Requirements

- R (recommended ≥ 4.1)
- Key packages commonly used in this project:
  - `limma`
  - `e1071` (Naive Bayes)
  - `pROC` (AUC)
  - `ggplot2`, `dplyr`, `tidyr` (plots + data handling)

Install packages (example):

```r
install.packages(c("e1071","pROC","ggplot2","dplyr","tidyr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("limma")
```

## How to run

> Scripts may be organized differently depending on the current repo structure. A typical workflow is:

1. **Run simulations** to reproduce FDR/power and AUC figures.
2. **Run real-data analysis** to produce the significant protein list and cross-validation AUC results.

If you want, paste your repo file tree (or tell me主要脚本文件名/路径), I can adapt this README to point to the *exact* script names and one-command run steps.

## Method summary (very short)

- Fit protein-wise linear models with an appropriate design matrix for group × time.
- Use empirical Bayes to stabilize protein-wise variance estimates.
- Test the global null hypothesis of no group difference across all time points using a moderated F-statistic.
- Adjust p-values using Benjamini–Hochberg (BH) to control FDR.
- Train Naive Bayes classifiers using top-ranked proteins; evaluate with subject-wise cross-validation.

## Citation

If you use this code, please cite the corresponding manuscript (update with DOI/journal info when available).

## License

Add a license for your preferred usage (e.g., MIT, GPL-3).  
Currently: not specified.
