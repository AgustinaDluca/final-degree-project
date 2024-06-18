**This folder contains all scripts with final analysis once T cells were removes in merged data containing all four samples. Here there is a brief explanation of the scripts:**

###differential abundance_testing.R: This script integrates employs sccomp_glm and miloR packages to explore differential abundance between two conditions (good prognosis and bad prognosis) and expression patterns in single-cell RNA sequencing data, focusing on identifying markers associated with the different conditions (prognosis).

###final_RNA_analysis.R: This script contains the final transcriptomics (RNA) analysis of the data, after the removal of T cells and all pertinent changes. The script systematically processes RNA sequencing data, from initial data loading and preprocessing to advanced analyses like clustering, integration, and marker gene identification, providing comprehensive insights into cellular composition and gene expression patterns for CLL prognosis.

###final_ATAC_analysis.R: This script contains the final epigenomics (RNA) analysis of the data, after the removal of T cells and all pertinent changes. The script systematically processes ATAC sequencing data, from initial data loading and preprocessing to advanced analyses like clustering, integration, and differentially accesible peaks identification, accesibility visualizations, providing comprehensive insights into chromatin accesibiloity patters and crucial peaks for CLL prognosis

###signature.R: This script is designated for the bad prognosis signature contruction, with the UCell package, and its visualizations. It focuses on exploring and visualizing gene expression patterns associated with specific gene sets related to bad prognosis in CLL in different clusters derived from the multiome analysis.
