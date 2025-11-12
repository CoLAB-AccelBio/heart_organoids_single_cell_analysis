# Single-cell analysis of heart organoids

🫀 Pipeline for analysing scRNA-seq data from human pluripotent stem cell-derived heart epicardium-myocardium organoids.

**Associated manuscript:** *Recreating coronary vascularization and sympathetic innervation of myocardium on a human pluripotent stem cell-derived heart assembloid* (submitted to _Cell Reports_).

## Overview

This repository contains code and documentation used for the analysis of scRNA-seq data generated from human induced pluripotent stem cells (hiPSC)-derived epicardium-myocardium organoid (EMO) models. The study compares **VEGFA-stimulated** (V-EMO) and **PDGFBB-stimulated** (P-EMO) organoids, profiling transcriptomic changes at the single-cell level.

## Data summary

Two organoid samples were analysed:

- **VEGFA** (V-EMO)
- **PDGFBB** (P-EMO)

Single-cell libraries were prepared using the *10x Genomics Chromium Next GEM Single Cell 3′ Kit v3.1* and sequenced on an *Illumina NextSeq 2000*.  
Raw data were processed with *Cell Ranger v9.0.1*, and downstream analyses were performed using ***[scStudio](https://www.biorxiv.org/content/10.1101/2025.04.17.649161v1)***  and custom [R scripts](https://github.com/CoLAB-AccelBio/heart_organoids_single_cell_analysis/tree/main/scripts).


## **Data availability**

- **Raw FASTQ** **files** are available at Sequence Read Archive (SRA) under the [PRJNA1356374](https://www.ncbi.nlm.nih.gov/sra/PRJNA1356374)
- **Processed data**, including count matrices, cell-level metadata, and annotated Seurat objects, are available at Gene Expression Omnibus (GEO) (pending accession number). Count matrices and cell-level metadata are also available in [data](https://github.com/CoLAB-AccelBio/heart_organoids_single_cell_analysis/tree/main/data).


## Data analysis workflow
### Pre-processing

- Quality control, filtering of low-quality cells
- Removal of doublets using `scDblFinder`
- Library size normalization using `scran`

### Dimensionality reduction and clustering analysis

- Feature selection and PCA using `Seurat`
- Dimensionality reduction via UMAP and t-SNE
- Clustering with shared nearest neighbor (SNN) modularity optimization
- Cluster stability assessment using `clustree`

### Differential expression analysis

- Performed using `Seurat::FindMarkers` with the **MAST** framework
- Marker gene identification and visualisation

### Trajectory inference

- Conducted with `TSCAN`
- Pseudotime analysis to identify lineage-specific transcriptional programs

### Visualisation

- Cell type annotation, marker expression plots
- PCA, t-SNE and trajectory plots using `ggplot2`

---

### Software and algorithms

| Software    | Version | Source                                                                             |
| ----------- | ------- | ---------------------------------------------------------------------------------- |
| Cell Ranger | 9.0.1   | [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger/latest)    |
| scStudio    | -       | [DiseaseTranscriptomicsLab](https://github.com/DiseaseTranscriptomicsLab/scStudio) |
| R           | 4.1.2   | [R Project](https://www.r-project.org)                                             |
| R           | 4.3.3   | [R Project](https://www.r-project.org)                                             |
| Seurat      | 4.1.1   | [Satija Lab](https://satijalab.org/seurat/)                                        |
| scDblFinder | 1.8.0   | [Bioconductor](https://bioconductor.org/packages/scDblFinder)                      |
| scran       | 1.22.1  | [Bioconductor](https://bioconductor.org/packages/scran)                            |
| scater      | 1.22.0  | [Bioconductor](https://bioconductor.org/packages/scater)                           |
| TSCAN       | 1.44.0  | [Bioconductor](https://bioconductor.org/packages/TSCAN)                            |
| clustree    | 0.5.1   | [CRAN](https://cran.r-project.org/web/packages/clustree/index.html)                |

