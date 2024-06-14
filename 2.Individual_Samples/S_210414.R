---
title: "patient S_210414"
author: "Agustina De Luca"
date: "2024-01-15"
output: html_document
---

## SAMPLE3

# What we did was to first pass all the raw data into the cell ranger ARC in order to get the data we can work with. Once this has finished in the cluster I 
# downloaded the folder named 'filtered_feature_bc_matrix' where all the three different matrix outputs are located
  
library(dplyr)
library(Seurat)
library(patchwork)
library(PCAtools)
library(devtools)
library(hdf5r)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biovizBase)
library(matchSCore2)

# Load the Seurat object back into R
S3 <- readRDS("S3_new.rds")

# Assuming your Seurat object is named 'your_seurat_object'
saveRDS(S3, file = "new_S3/S3.rds")

### How to load multiome data

counts <- Read10X_h5("/home/adeluca/Documents/S_210414/filtered_feature_bc_matrix.h5")
names(counts)

data <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# Bash
# tabix -p bed /home/adeluca/Documents/S_210414/atac_fragments.tsv.gz

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#seqlevelsStyle(annotation) <- "UCSC"
ucsc.levels <- stringr::str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels
genome(annotation) <- "hg38"

data[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = "/home/adeluca/Documents/S_210414/atac_fragments.tsv.gz",
  annotation = annotation
)

## I checked that the the ZAP70 was present in the 'annotation' genome we selected (hg38)
## https://www.ncbi.nlm.nih.gov/datasets/gene/id/dc117a1eacd2c3839fde463efbb051fc/

dim(data)
# 36601  3390
head(rownames(data), 5) # genes
head(colnames(data), 5) # cells

data[["RNA"]]$counts

head(data[[]], 3)
# [1] "orig.ident"   "nCount_RNA"   "nFeature_RNA"
# [4] "percent.mt"  

## 1. QC-filtering of cells

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-") # percentage of mitochondrial genes in each cell

head(data@meta.data, 5)
### SAVE ###
# Visualize the different metrics with violin plots
VlnPlot(S3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

hist(S3$percent.mt, main = 'Distribution of mitochondrial genes in S3')
hist(S3$nFeature_RNA, main = 'Distribution of the nFeatures in S3')

summary(data[['nFeature_RNA']])
summary(data[['percent.mt']])

# Visualize feature-feature relationship with scatterplots
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
# We expect this plot to not have any type of correlation as plotting the depth of the samples with the percentage of mitochondrial RNA. In fact we wan more RNA without mitochondrial. The correlation coefficient is -0.08.
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# We expect this plot to have a strong correlation as we are relating the number of RNA counts (depth) and the number of features meaning the genes. This can be explained as the higher the depth the better the coverage and sensitivity in detecting gene expression making the analysis much more comprehensive analysis, moreover higher numbers of features suggest greater biological complexity and diversity in terms of the genes being expressed.

# Deeper sequencing will result in the detection of a greater number of expressed genes, reflecting the biological complexity of the cells being analyzed.

plot1+plot2

hist(data@meta.data$nFeature_RNA)

### SELECT THOSE CELLS THAT HAVE A GOOD SCORE REGARDING NUMBER OF FEATURES AND THE PERCENTAGE OF MITOCHONDRIAL CONTAMINATION!!!

#data <- subset(data, subset = nFeature_RNA > 500 & percent.mt < 10) # We are going to select the cells we want to be analysed based on the previous results mainly the violin plots from above
data <- subset(data, subset = nFeature_RNA > 500 & percent.mt < 20) # We are going to select the cells we want to be analysed based on the previous results mainly the violin plots from above

# I could put a threshold in the percent.mt in 5 however it is not needed as the values as low it doesn't overcome the more than the 20 percent of the samples

### SAVE ###
# Visualize the different metrics with violin plots
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

# Visualize feature-feature relationship with scatterplots
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

summary(data@meta.data$nFeature_RNA) # It is really interesting to look at the summary of the number of features as it helps us to choose a threshold to remove problematic cells
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 504.0   898.2  1359.0  1717.0  1901.8  9806.0

## With the filter in < 20 mitochondrial genes
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 502     871    1279    1644    1825    9806 

summary(data@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3785  5.9280  7.4966  7.2069  8.8180  9.9929 

## With the filter in < 20 mitochondrial genes
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3785  7.4877  9.9456 10.4849 13.4076 19.9664 
#plot(density(log10(data@meta.data$nFeature_RNA)))

## 2. NORMALIZATION

data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

data[["RNA"]]$counts

## 3. IDENTIFY HIGHL VARIABLE FEATURES (feature selection)

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000) # 15-20% of the number of genes  TRY WITH 800

top10 <- head(VariableFeatures(data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## 4. SCALING - ensures that the magnitudes of gene expression are comparable across all genes, preventing genes with higher magnitudes from dominating the analysis.

data <- ScaleData(data) # This function scales and centers the expression values for each gene across all cells. 

#head(data[['RNA']]$scale.data, 5)

## 5. LINEAR DIMENSIONALITY REDUCTION (PCA)

data <- RunPCA(data, features = VariableFeatures(object = data))

print(data[['pca']], dims = 1:5, nfeatures = 5) # Examine and visualize PRINCIPAL COMPONENT results a few different ways

VizDimLoadings(data, dims = 1:2, reduction = 'pca') # shows the results of the PC's and the genes that are positive and negative. This function visualizes the loadings of features on the first two principal components. It helps you understand which genes contribute the most to the observed variation along these components.

DimPlot(data, reduction = 'pca') + NoLegend() # The position of cells along PC1 and PC2 represents their overall gene expression patterns with respect to the main sources of variation in the dataset. 
### SAVE ###
DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)
### SAVE ###
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)

## 6. DETERMINE THE DIMENSIONALITY OF THE DATASET

ElbowPlot(data)

'''
The variable stdev contains a vector of the PCs standard deviations
#explained_variance <- (data[["pca"]]@stdev)^2
#barplot(explained_variance)

#percentage <- (explained_variance / sum(explained_variance)) * 100

#BiocManager::install("PCAtools")
#library(PCAtools)

#optimal_PC <- findElbowPoint(percentage) # find the optimal PC to retain for the subsequent analyses by looking for the elbow point in the plot

#plot(percentage, xlab="PC", ylab="Variance explained (%)")
#abline(v=optimal_PC, col="red")
'''

## 7. CLUSTER THE CELLS

data <- FindNeighbors(data, dims = 1:20) #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.

data <- FindClusters(data, resolution = 0.5) #Partition the KNN graph into clusters.

table(data@active.ident)
# 0   1   2   3   4   5 
# 371 217 201 166 103  16 

## With the filter < 20 in mitochondrial genes
# 0   1   2   3   4   5   6 
# 545 465 313 304 270 137  95 
head(Idents(data), 5) #Cluster assignments for each cell after the clustering process.

## 8. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

data <- RunUMAP(data, dims = 1:20)
### SAVE ###
DimPlot(S3, reduction = 'umap', label = T, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','lightskyblue', 'violet'))

FeaturePlot(data, features = 'nCount_RNA')
FeaturePlot(data, features = 'nFeature_RNA')
FeaturePlot(data, features = 'percent.mt')

### SAVE ###
VlnPlot(data,features = "nFeature_RNA",log = T, sort = T, c('violet','#F4D35E','#118AB2','#440154ff','#FF6B6B','#8338EC', '#06D6A0'))

## 9. FINDING DIFFERENTIALLY EXPRESSED FEATURES (cluster biomarkers)

markers <- FindAllMarkers(S3, only.pos = T)
### SAVE ###
# Save my markers matrix
#write.csv(markers, file = "new_S3/markers_S3_new.csv", row.names = T)
markers <- read.csv('/home/adeluca/Documents/markers_S3_new.csv', header = F, col.names = T)
names(markers) <- c('markers','p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'cluster', 'gene')

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(data, features = top10$gene) + NoLegend()


library(matchSCore2)
top <- cut_markers(clusters = levels(as.factor(markers$cluster)),markers, ntop=50)
top
### SAVE ###
#write.csv(top, file = "new_S3/top50markers_S3_new.csv", row.names = T)

### SAVE ###
#Most common markers
FeaturePlot(data, features = c("MS4A1","GNLY", "CD3E", "CD14", "FCE51A", "FCGR3A", "LYZ", "PPBP", "CD8A"),order = T) # location in the clusters
VlnPlot(data, features = c("MS4A1","GNLY", "CD3E", "CD14", "FCE51A", "FCGR3A", "LYZ", "PPBP", "CD8A"), slot = "counts", log = TRUE) # RNA counts

mark3 <- read.csv('new_S3/top50markers_S3_new.csv')

#Cluster0
FeaturePlot(data, features = c("CHRNA10", "AC073934.1", "KIRREL2", "AC005329.2", "DNAJB5-DT", "RRH", "AC083799.1", "OSGEPL1-AS1", "AC096759.1", "AC010983.1"), order = T) # location in the clusters
FeaturePlot(data, features = c("CHRNA10", "KIRREL2", "AC005329.2", "RRH"), order = T) # location in the clusters
FeaturePlot(data, features = "CHRNA10", order = T) # location in the clusters
FeaturePlot(data, features = "KIRREL2", order = T) # location in the clusters
FeaturePlot(data, features = "AC005329.2", order = T) # location in the clusters
FeaturePlot(data, features = "RRH", order = T) # location in the clusters


FeaturePlot(data, features = c("LINC01361", "AL358933.1", "AC011468.1", "AL596218.1", "DNAJB5-DT", "GLMP", "LINC01953", "AC116099.1", "AL133230.1", "AC010983.1"), order = T) # location in the clusters

#Cluster1
FeaturePlot(data, features = c("GALNT16", "AC007610.2", "LINC02774", "AC095057.3", "AL022329.1", "AC007533.1", "GUCY2C", "SHC3", "HIST1H2BB", "FHL1"), order = T) # location in the clusters
FeaturePlot(data, features = c("GALNT16", "AC007610.2", "HIST1H2BB", "FHL1"), order = T) # location in the clusters
FeaturePlot(data, features = "GALNT16", order = T) # location in the clusters
FeaturePlot(data, features = "AC007610.2", order = T) # location in the clusters
FeaturePlot(data, features = "HIST1H2BB", order = T) # location in the clusters
FeaturePlot(data, features = "FHL1", order = T) # location in the clusters

FeaturePlot(data, features = c("IGKC", "PPOX", "HLA-DPA1", "IGHM", "MMP24", "RCN3", "BCAT2", "CSMD1", "GHDC", "ZBTB34"), order = T) # location in the clusters

#Cluster2
FeaturePlot(data, features = c("PLA2R1", "CARMIL3", "PRDM12", "ATP11A-AS1", "SHROOM4", "C8orf86", "SEC62-AS1", "IZUMO4", "VPS28", "IGKC"), order = T) # location in the clusters
FeaturePlot(data, features = c("PLA2R1", "CARMIL3", "PRDM12", "IGKC"), order = T) # location in the clusters
FeaturePlot(data, features = "PLA2R1", order = T) # location in the clusters
FeaturePlot(data, features = "CARMIL3", order = T) # location in the clusters
FeaturePlot(data, features = "PRDM12", order = T) # location in the clusters
FeaturePlot(data, features = "IGKC", order = T) # location in the clusters

FeaturePlot(data, features = c("AC017104.6", "RTP5", "KLRC4", "TNFRSF25", "PTGDR", "TNFRSF18", "LINC02084", "LINC01943", "KIF19", "AC021028.1"), order = T) # location in the clusters

#Cluster3
FeaturePlot(data, features = c("CACNA1C", "TRPM3", "GALNT17", "KCNMB2", "AC073114.1", "NTM", "OBI1-AS1", "LINC01122", "RP1", "LINC02267"), order = T) # location in the clusters
FeaturePlot(data, features = c("CACNA1C", "GALNT17", "KCNMB2", "LINC02267"), order = T) # location in the clusters
FeaturePlot(data, features = "CACNA1C", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "GALNT17", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "KCNMB2", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "LINC02267", order = T, min.cutoff = 1) # location in the clusters

FeaturePlot(data, features = c("FAIM2", "SLC17A1", "AL132657.1", "RLBP1", "AC087393.2", "LRRC73", "AL513128.3", "AC100835.2", "RNF208", "AC011472.4"), order = T) # location in the clusters

#Cluster4
FeaturePlot(data, features = c("TNFRSF25", "TNFRSF18", "AC017104.6", "LINC01973", "LGALS3BP", "BEX3", "XCL1", "ULBP2", "WNT10B", "KRT72"), order = T) # location in the clusters
FeaturePlot(data, features = c("TNFRSF25", "TNFRSF18", "XCL1", "ULBP2"), order = T) # location in the clusters
FeaturePlot(data, features = "TNFRSF25", order = T) # location in the clusters
FeaturePlot(data, features = "TNFRSF18", order = T) # location in the clusters
FeaturePlot(data, features = "XCL1", order = T) # location in the clusters
FeaturePlot(data, features = "ULBP2", order = T) # location in the clusters

FeaturePlot(data, features = c("AC007610.2", "AL022329.1", "TMEM176A", "EFNB2", "FHL1", "HLX-AS1", "PLD5", "AC005699.1", "AL117339.4", "L3MBTL4-AS1"), order = T) # location in the clusters

#Cluster5
FeaturePlot(data, features = c("CD300LB", "SIRPA", "P2RY2", "VENTX", "TNFAIP2", "AL355881.1", "MILR1", "AL354928.1", "AL356056.2", "AC125603.1", "MS4A7"), order = T) # location in the clusters
FeaturePlot(data, features = c("CD300LB", "SIRPA", "MS4A7", "TNFAIP2"), order = T) # location in the clusters
FeaturePlot(data, features = "CD300LB", order = T) # location in the clusters
FeaturePlot(data, features = "SIRPA", order = T) # location in the clusters
FeaturePlot(data, features = "MS4A7", order = T) # location in the clusters
FeaturePlot(data, features = "TNFAIP2", order = T) # location in the clusters

FeaturePlot(data, features = c("CLEC4E", "AC104809.2", "HES1", "ZNF697", "AC093627.3", "P2RY2", "BTBD3", "ASGR2", "TNFAIP2", "CLEC11A"), order = T) # location in the clusters

#Cluster6

FeaturePlot(data, features = c("AC119868.2", "RIT2", "ZNF536", "AC008571.2", "AC004946.1", "RGS7BP", "AL365214.2", "PAK5", "LINC00348", "FIGN"), order = T) # location in the clusters

# Color each cluster

DimPlot(data, group.by = "seurat_clusters", cols = c('#118AB2','grey', 'grey', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','#06D6A0', 'grey', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', '#8338EC', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', '#F4D35E', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', '#FF6B6B', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'grey', '#440154ff', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'grey', 'grey', 'violet'))


grep('IL7R', mark$cluster6)
grep('CCR7', mark$cluster6)
grep('S100A4', mark$cluster6)
grep('CD14', mark$cluster6)
grep('LYZ', mark$cluster6)
grep('MS4A1', mark$cluster6)
grep('CD8A', mark$cluster6)
grep('FCGR3A', mark$cluster6)
grep('MS4A7', mark$cluster6)
grep('GNLY', mark$cluster6)
grep('NKG7', mark$cluster6)
grep('PPBP', mark$cluster6)
grep('FCER1A', mark$cluster6)
grep('CST3', mark$cluster6)
grep('HBA', mark$cluster6)
grep('HBB', mark$cluster6)

# Subclustering

Idents(S3) <- S3$seurat_clusters
S3 <- FindNeighbors(S3, dims = 1:20, graph.name = "graph") #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.
S3 <- FindSubCluster(S3, resolution = 0.2, cluster = '2', graph.name = "graph")
table(S3$sub.cluster)

Idents(S3) <- as.factor(S3$sub.cluster)

DimPlot(S3, reduction = "umap", label = T, pt.size = 0.5,  cols = c('#118AB2','#06D6A0','#8338EC','orange','#F4D35E', '#FF6B6B', 'lightskyblue', 'violet')) + NoLegend()

markers3 <- FindAllMarkers(S3, only.pos = T)
top3 <- cut_markers(clusters = levels(markers3$cluster),markers3, ntop=50)

write.csv(markers3, 'new_S3/markers_S3_new.csv')
write.csv(top3, 'new_S3/top_markers_subclustering_S3.csv')

### CELLTYPE ANNOTATION (AFTER SEARCHING FOR SPECIFIC MARKERS IN CELLTYPIST)

S3 <- SetIdent(S3, value = as.factor(S3$sub.cluster))

cluster_annotation <- c("Leukemic cells 1", "Leukemic cells 2", "CD4+ T cells", "CD8+ Cytotoxic T cells", "Leukemic cells 3 (IG+)", "Leukemic cells 4", "Myeloid cells", "Leukemic cells 5")
cell_type <- cluster_annotation[as.factor(levels(S3))]
S3$cell_type <- cell_type[as.factor(S3$sub.cluster)]

names(cluster_annotation) <- levels(S3)

S3 <- SetIdent(S3, value = as.factor(S3$cell_type))

DimPlot(S3, reduction = "umap", label = T, pt.size = 0.5, cols = c('violet','pink','mediumseagreen','mediumaquamarine','darkseagreen','seagreen', 'green4', '#F4D35E')) + NoLegend()
DimPlot(S3, reduction = "umap", label = F, pt.size = 0.5, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','#440154ff')) + NoLegend()

## SHINY APP
setwd("/home/adeluca/Documents/")
scConf = createConfig(S3)
makeShinyApp(S3, scConf, gene.mapping = TRUE,
             shiny.title = "Sample3")

# SIMULACION - Azimuth Seurat











