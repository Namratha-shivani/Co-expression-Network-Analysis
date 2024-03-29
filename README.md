# Co-expression Network Analysis
Co-expression network analysis is a computational approach used in bioinformatics and systems biology to study the relationships between genes based on their expression patterns across multiple samples or conditions. It is helpful for

1. Identifying Functional Relationships
2. Understanding Regulatory Mechanisms
3. Biomarker Analysis
4. Comparitive Analysis
5. Integration with omics data

# Project 

<h1>_Exploration of Common and Differentially Expressed Genes in Psoriasis, and Dermatomyositis_</h1>

### Introduction 

Psoriasis, and dermatomyositis (DM) are autoimmune inflammatory diseases that affect millions of people worldwide. Psoriasis is characterized by skin manifestations, whereas DM affects skeletal muscles and skin. Despite their differences, these diseases share common immune dysregulation pathways and treatment modalities.

### Objective 

This project aims to identify common and differentially expressed genes (DEGs) among psoriasis, and DM using weighted gene coexpression network analysis (WGCNA). By exploring gene expression patterns, we seek to uncover potential treatment targets and biomarkers shared among these diseases.

# Weighted Gene Co-expression Network Analysis

The WGCNA R software package is a toolbox of tools for analyzing gene expression data. It helps researchers understand how genes work together by looking at patterns in their activity levels. With WGCNA, you can build networks to see which genes are connected, find groups of genes that behave similarly, and visualize your results to make sense of complex data.

```r
# Install Packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install('WGCNA')

library(WGCNA)

```
