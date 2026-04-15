# Harmonisation-Paper

This repository contains code and data for running, evaluating, and comparing harmonisation methods designed to remove batch effects from T1-weighted (T1w) brain MRI–derived outcome measures (imaging-derived phenotypes; IDPs).

The work focuses on harmonisation across site, scanner, or mixed batch effects using the DPUK PET-MR test–retest dataset.

---

## Overview

Multi-site neuroimaging studies are susceptible to systematic variability introduced by differences in scanners, acquisition protocols, and imaging sites. This project provides a framework for:

- Apply multiple harmonisation approaches  
- Quantitatively evaluate harmonisation performance  
- Compare the impact of harmonisation on imaging-derived phenotypes (IDPs)  

The aim is to support transparent, reproducible, and robust assessment of batch-effect correction methods.

---

## Repository Structure

```text
Harmonisation-Paper/
│
├── harmonisation_methods/
│   Codes used to run harmonisation methods on the original (raw) IDPs
│
├── harmonisation_diagnostics/
│   Codes and functions used to calculate harmonisation evaluation metrics
│
├── manuscript materials/
│   Appendix A: IQM-informed harmonisation methodology
│   supplementary docx and excel file with additional results
│
└── README.md

```

## Harmonisation Methods

The `harmonisation_methods/` directory includes scripts used to:

- Apply harmonisation techniques to raw IDPs  
- Remove batch effects associated with site, scanner, or mixed sources  

These implementations can be adapted for use in other multi-site neuroimaging studies.

---

## Harmonisation Diagnostics

The `harmonisation_diagnostics/` directory contains code to:

- Quantify residual batch effects after harmonisation  
- Evaluate preservation of biological signal  
- Compare performance across harmonisation methods  

All metrics are designed to support objective and reproducible evaluation.

---

## Use Cases

This repository is intended for:

- Methodological research on batch-effect correction  
- Benchmarking harmonisation approaches  
- Reproducible analysis of multi-site neuroimaging data  

## Python-version

- The evaluation metrics in this study are now included in the DiagnoseHarmonise library. It is written fully in Python (See [DHARM](https://jake-turnbull.github.io/HarmonisationDiagnostics/)) 
---

## Contact

For questions, feedback, or contributions, please open an issue or contact [Gaurav Bhalerao](gaurav.bhalerao@ndcn.ox.ac.uk).
