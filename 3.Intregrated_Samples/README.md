**This folder contains all analysis related to merged samples with the T cells. 
Here you will find a brief description of each code:**

**RNA_ATAC_all_samples.R**: This script outlines a robust workflow for analyzing multi-omics single-cell data, specifically focusing on integrating RNA and ATAC-seq profiles. It employs Seurat in R to conduct a comprehensive analysis pipeline that includes data loading, quality control (QC), normalization, dimensionality reduction (like PCA and t-SNE), clustering of cells based on similarities in their ATAC-seq profiles, differential accessibility analysis to identify significant chromatin changes across conditions or cell types, and motif enrichment analysis to uncover enriched DNA-binding motifs.

**RNA_all_samples.R**: This script is designed for analyzing single-cell RNA-seq data from all merged samples using Seurat in R. It follows a structured pipeline that includes data loading, quality control (QC), normalization, dimensionality reduction (such as PCA and t-SNE), cell clustering based on gene expression profiles, identification of marker genes for different cell types or states, and visualization of results. This approach helps uncover cellular heterogeneity and biological insights from RNA-seq data, enabling detailed exploration of gene expression patterns across individual cells or conditions.


