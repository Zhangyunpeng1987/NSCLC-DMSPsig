# NSCLC-DMSP.sig

## Overview

This repository provides the R code used in the study:

"Coordinated multicellular immune programs and drug targets revealed by
single-cell analysis in driver-mutated NSCLC"

The project aims to systematically characterize tumor immune
microenvironment (TIME) heterogeneity in driver-mutated non-small cell
lung cancer (NSCLC) and to construct a prognostic model (DMSP.sig).

------------------------------------------------------------------------

## Main Analyses Included

-   Single-cell RNA-seq integration (Seurat + harmony)
-   Cell type annotation and subclustering
-   TIME module identification (CM1-CM5)
-   Differential expression and enrichment analysis
-   Gene set scoring
-   Cancer cell state analysis (copykat + infercnv)
-   TF activity (decoupleR)
-   Cell-cell communication (CellChat)
-   Metabolic activity analysis (scMetabolism)
-   Survival analysis (TCGA)
-   Prognostic model construction (LASSO Cox)
-   Immune infiltration analysis
-   Drug sensitivity prediction

------------------------------------------------------------------------

## Code Organization

The scripts are organized by figure to facilitate reproducibility. Each folder contains the code used to generate the corresponding panels in the manuscript.

For example:
- `code/Figure1/` contains scripts for Figure 1
- `code/FigureS1/` contains scripts for supplementary Figure 1.
- Subsequent folders correspond to the remaining figures.

Users can follow these scripts to reproduce the results step by step.

------------------------------------------------------------------------

## Repository Structure

NSCLC-DMSPsig/
├── README.md
├── LICENSE
├── code/
│   ├── Figure1/
│   │   ├── Figure1.R
│   │
│   ├── Figure2/
│   │   ├── Figure2.R
│   │
│   ├── Figure3/
│   │   ├── Figure3.R
│   │
│   ├── Figure4/
│   │   ├── Figure4.R
│   │
│   ├── Figure5/
│   │   ├── Figure5.R
│   │
│   ├── Figure7/
│   │   ├── Figure7.R
│   │
│   ├── FigureS1/
│   │   ├── FigureS1.R
│   │
│   ├── FigureS2/
│   │   ├── FigureS2.R
│   │
│   ├── FigureS3/
│   │   ├── FigureS3.R
│   │
│   ├── FigureS4/
│   │   ├── FigureS4.R
│   │
│   ├── FigureS5/
│   │   ├── FigureS5.R
│   │
│   ├── FigureS6/
│   │   ├── FigureS6.R
│   │
│   ├── FigureS7/
│   │   ├── FigureS7.R
│
└── sessionInfo.txt

------------------------------------------------------------------------

## Data Availability

The datasets used in this study are publicly available:

- Single-cell RNA-seq datasets were obtained from the Gene Expression Omnibus (GEO) database (https://www.ncbi.nlm.nih.gov/geo/), with accession numbers provided in the manuscript.
- TCGA bulk RNA-seq data: https://portal.gdc.cancer.gov/
- Spatial transcriptomics data: E-MTAB-13530 (EMBL-EBI)
- Human Protein Atlas: https://www.proteinatlas.org/

Due to file size limitations, raw datasets are not included in this repository.

------------------------------------------------------------------------

## Requirements

R = 4.4.2

Packages: Seurat, harmony, GSVA, glmnet, CellChat, scMetabolism, decoupleR, copykat,
infercnv, oncoPredict

------------------------------------------------------------------------

## Contact

Yunpeng Zhang: zhangyp@hrbmu.edu.cn\
Xia Li: lixia@hrbmu.edu.cn

------------------------------------------------------------------------

## License

MIT License
