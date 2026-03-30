# 🧬 Transcriptomic Analysis of Breast Cancer Using TCGA RNA-seq (DESeq2)

## Project Goal
Developing a transcriptomic analysis pipeline to perform differential gene expression analysis on TCGA breast cancer RNA-seq data and identify subtype-associated genes and biological pathways.

---

## Overview
In this project, I am building a bulk RNA-seq analysis workflow using TCGA breast cancer data to identify differentially expressed genes (DEGs) and characterize subtype-specific transcriptional patterns.

This pipeline reflects standard transcriptomic workflows used in cancer genomics research.

---

## Workflow Overview

1. Data preprocessing and normalization  
2. Differential gene expression analysis (DESeq2)  
3. Visualization (PCA, volcano plot)  
4. Functional enrichment analysis (GO / pathways)  

---

## Key Features

- Differential gene expression analysis using DESeq2  
- Identification of subtype-associated genes  
- PCA-based visualization of sample variation  
- Functional enrichment analysis using clusterProfiler  
- Reproducible and modular analysis pipeline  

---

## Project Structure

tcga-rna-seq-breast-cancer/  
│── scripts/  
│── data/  
│── results/  
│── figures/  
│── README.md  

---

## Workflow Details

### 1️⃣ Data Preprocessing
- Load TCGA RNA-seq count matrix and clinical metadata  
- Normalize gene expression counts  

### 2️⃣ Differential Expression Analysis
- Perform DESeq2 analysis across breast cancer subtypes  
- Identify significantly differentially expressed genes  

### 3️⃣ Visualization
- PCA plot for sample clustering  
- Volcano plot for differential expression  

### 4️⃣ Functional Enrichment
- Gene Ontology (GO) enrichment  
- Pathway analysis using clusterProfiler  

---

## Status

This project is currently in progress. The analysis pipeline structure is implemented, and differential expression and enrichment analyses are being developed.

---

## Expected Outputs

- List of differentially expressed genes  
- PCA plot showing subtype separation  
- Volcano plot of gene expression changes  
- Enriched biological pathways  

---

## Tools & Technologies

- R: DESeq2, clusterProfiler, ggplot2  
- Python: pandas, NumPy  

---

## Skills Demonstrated

- Bulk RNA-seq analysis  
- Differential gene expression  
- Cancer genomics analysis  
- Functional enrichment analysis  
- Data visualization  
- Bioinformatics pipeline development  

---

## Impact

This project demonstrates transcriptomic analysis of breast cancer data and reflects workflows used in:

- Cancer genomics research  
- Biomarker discovery  
- Precision medicine  

---

## Author

Divya Reddy  
MS Bioinformatics, Georgia Institute of Technology  
