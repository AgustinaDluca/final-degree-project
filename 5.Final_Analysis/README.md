This folder contains all scripts with final analysis once T cells were removes in merged data.

differential abundance_testing.R: This script integrates employs sccomp_glm and miloR packages to explore differential abundance between two conditions (good prognosis and bad prognosis) and expression patterns in single-cell RNA sequencing data, focusing on identifying markers associated with the different conditions (prognosis).

final_RNA_analysis.R: Final transcriptomics analysis with all merged samples (RNA Assay)
final_RNA_analysis.R: Final epigenomics analysis with all merged samples (ATAC Assay)

signature.R:  this script is designated for the bad prognosis signature contruction, with the UCell package, and its visualizations. It focuses on exploring and visualizing gene expression patterns associated with specific gene sets related to bad prognosis in CLL in different clusters derived from the multiome analysis.
