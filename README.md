# scRNA-seq Analysis for MSc Research Project

## Overview

This repository contains the code for the scRNA-seq analysis conducted as part of my MSc research project. The analysis is divided into two main pipelines:

1. **Upstream Analysis Pipeline**
2. **Downstream Analysis Pipeline**

These pipelines are designed to ensure efficient and reproducible workflows for processing and analysing single-cell RNA sequencing (scRNA-seq) data.

## 1. Upstream Analysis Pipeline

The upstream analysis pipeline focuses on the initial processing and preparation of the raw scRNA-seq data. The steps included in this pipeline are:

1. **Raw Data (Count Table):**  
   The pipeline begins with the raw count table obtained from sequencing.
   
2. **Data Cleaning:**  
   This step involves filtering and cleaning the raw data to remove any irrelevant or low-quality entries.
   
3. **Quality Control:**  
   Quality control checks are performed to ensure the data is of high quality, removing any cells or genes that do not meet specific criteria.
   
4. **Doublet Removal:**  
   This step identifies and removes potential doublets, which are instances where two cells are captured together, to avoid skewing the analysis.
   
5. **Normalisation:**  
   The data is normalised to ensure that differences in sequencing depth or library size do not affect downstream analysis.
   
6. **Feature Selection:**  
   Relevant features (genes) are selected based on variability across cells to be used in further analysis.
   
7. **PCA (Principal Component Analysis):**  
   PCA is performed to reduce the dimensionality of the data, helping to identify major trends and patterns.
   
8. **kNN/UMAP:**  
   Nearest neighbor clustering (kNN) and UMAP are applied to visualise the data in lower dimensions.
   
9. **Clustering:**  
   Cells are clustered based on similarity, identifying groups of cells with similar expression profiles.
   
10. **Cell Annotation:**  
    Finally, the clusters are annotated to assign biological meaning, identifying different cell types or states.

The output of this pipeline is an annotated data table that serves as the input for downstream analysis.

## 2. Downstream Analysis Pipeline

The downstream analysis pipeline focuses on deriving biological insights from the annotated data table produced by the upstream pipeline. The steps include:

1. **DEG (Differential Expression Gene) Analysis:**  
   Identification of genes that are differentially expressed between conditions or clusters.
   
2. **Gene Ontology (GO) Analysis:**  
   GO analysis is performed to identify biological processes, cellular components, and molecular functions that are enriched in the differentially expressed genes.
   
   - **GO Enrichment Map:**  
     Visualisation of the GO terms enriched in the dataset.
   
3. **KEGG Enrichment Analysis:**  
   The KEGG pathway analysis identifies pathways that are enriched in the differentially expressed genes.
   
   - **Manual KEGG Pathway Analysis:**  
     Further manual curation and interpretation of the KEGG pathways to understand the underlying biological processes.

This pipeline provides insights into the biological functions and pathways associated with the identified cell populations.

---

## How to Use

To run the analysis, follow the instructions in the respective Jupyter notebooks and R scripts available in the repository. The upstream pipeline is implemented in Python using the Scanpy library, while the downstream pipeline is executed in R, leveraging tools such as Limma for DEG analysis and other R packages for pathway enrichment.

---
