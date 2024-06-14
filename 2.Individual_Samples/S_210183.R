---
title: "patient S_210183"
author: "Agustina De Luca"
date: "2024-01-15"
output: html_document
---
## SAMPLE2

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
S2 <- readRDS("S2.rds")

# Assuming your Seurat object is named 'your_seurat_object'
saveRDS(S2, file = "S2.rds")

### How to load multiome data

counts <- Read10X_h5("../adeluca/Documents/Documents/scRNA_samples/S_210183/filtered_feature_bc_matrix.h5")
names(counts)

data2 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# Bash
# tabix -p bed /home/adeluca/Documents/S_210183/atac_fragments.tsv.gz

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#seqlevelsStyle(annotation) <- "UCSC"
ucsc.levels <- stringr::str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels
genome(annotation) <- "hg38"

data2[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = "../adeluca/Documents/Documents/scRNA_samples/S_210183/atac_fragments.tsv.gz",
  annotation = annotation
)

## I checked that the the ZAP70 was present in the 'annotation' genome we selected (hg38)
## https://www.ncbi.nlm.nih.gov/datasets/gene/id/dc117a1eacd2c3839fde463efbb051fc/

dim(data)
# 36601  8084
head(rownames(data), 5) # genes
head(colnames(data), 5) # cells

data[["RNA"]]$counts

head(data[[]], 3)
# [1] "orig.ident"   "nCount_RNA"   "nFeature_RNA"
# [4] "percent.mt"  

## 1. QC-filtering of cells

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-") # percentage of mitochondrial genes in each cell

head(data@meta.data, 5)

# Visualize the different metrics with violin plots
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

hist(S2$percent.mt, main = 'Distribution of mitochondrial genes in S2')
hist(S2$nFeature_RNA, main = 'Distribution of the nFeatures in S2')

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

data <- subset(data, subset = nFeature_RNA > 500 & percent.mt < 10) # We are going to select the cells we want to be analysed based on the previous results mainly the violin plots from above

# I could put a threshold in the percent.mt in 5 however it is not needed as the values as low it doesn't overcome the more than the 20 percent of the samples

# Visualize the different metrics with violin plots
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

# Visualize feature-feature relationship with scatterplots
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

summary(data@meta.data$nFeature_RNA) # It is really interesting to look at the summary of the number of features as it helps us to choose a threshold to remove problematic cells
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 505    1108    1294    1647    1676    9783 

summary(data@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1668  2.2407  3.0729  3.4576  4.2673  9.9765 

#plot(density(log10(data@meta.data$nFeature_RNA)))

## 2. NORMALIZATION

data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

data[["RNA"]]$counts

## 3. IDENTIFY HIGHL VARIABLE FEATURES (feature selection)

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000) # 15-20% of the number of genes 

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

DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)

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
#   0    1    2    3    4    5    6 
# 2264 2021 1451 1281  456  322   25 

head(Idents(data), 5) #Cluster assignments for each cell after the clustering process.

## 8. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

data <- RunUMAP(data, dims = 1:20)

DimPlot(data, reduction = 'umap', label = T, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','#440154ff', 'violet'))

FeaturePlot(data, features = 'nCount_RNA')
FeaturePlot(data, features = 'nFeature_RNA')
FeaturePlot(data, features = 'percent.mt')

VlnPlot(data,features = "nFeature_RNA",log = T,sort = T)

## 9. FINDING DIFFERENTIALLY EXPRESSED FEATURES (cluster biomarkers)

markers <- FindAllMarkers(data, only.pos = T)

# Save my markers matrix
#write.csv(markers, file = "/home/adeluca/Documents/S_210183/markers_S210183.csv", row.names = T)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(data, features = top10$gene) + NoLegend()

library(matchSCore2)
top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=50)
top

#write.csv(top, file = "/home/adeluca/Documents/S_210183/top50markers_S210183.csv", row.names = T)

#Most common markers
FeaturePlot(data, features = c("MS4A1","GNLY", "CD3E", "CD14", "FCE51A", "FCGR3A", "LYZ", "PPBP", "CD8A"),order = T) # location in the clusters
VlnPlot(data, features = c("MS4A1","GNLY", "CD3E", "CD14", "FCE51A", "FCGR3A", "LYZ", "PPBP", "CD8A"), slot = "counts", log = TRUE) # RNA counts

mark2 <- read.csv('/home/adeluca/Documents/S_210183/markers/top_markers_S210183.csv')

#Cluster0
FeaturePlot(data, features = c("LINC01239", "AP003352.1", "PIPOX", "ZKSCAN5", "AC068051.1", "MRPL45", "AC022509.2", "CD24", "ZNF460-AS1", "SCART1"), order = T) # location in the clusters
FeaturePlot(data, features = c("LINC01239", "AP003352.1", "CD24", "SCART1"), order = T) # location in the clusters

FeaturePlot(data, features = "LINC01239", order = T) # location in the clusters
FeaturePlot(data, features = "AP003352.1", order = T) # location in the clusters
FeaturePlot(data, features = "CD24", order = T) # location in the clusters
FeaturePlot(data, features = "SCART1", order = T) # location in the clusters

#Cluster1
FeaturePlot(data, features = c("STK32B", "SGSM1", "WDR17", "P2RY14", "CDKL1", "MYOM1", "C1orf159", "KLHL14", "AC090125.1", "AUTS2"), order = T) # location in the clusters
FeaturePlot(data, features = c("STK32B", "SGSM1", "WDR17", "P2RY14"), order = T) # location in the clusters

FeaturePlot(data, features = "STK32B", order = T) # location in the clusters
FeaturePlot(data, features = "SGSM1", order = T) # location in the clusters
FeaturePlot(data, features = "WDR17", order = T) # location in the clusters
FeaturePlot(data, features = "P2RY14", order = T) # location in the clusters

#Cluster2
FeaturePlot(data, features = c("AC087627.1", "AC112484.1", "POF1B", "WNT6", "AC005586.1", "AC013287.1", "AC034102.6", "AL031710.1", "AL136980.1", "AC005699.1"), order = T) # location in the clusters
FeaturePlot(data, features = c("AC087627.1", "AC112484.1", "AC005586.1", "AC013287.1"), order = T) # location in the clusters

FeaturePlot(data, features = "AC087627.1", order = T) # location in the clusters
FeaturePlot(data, features = "AC112484.1", order = T) # location in the clusters
FeaturePlot(data, features = "AC005586.1", order = T) # location in the clusters
FeaturePlot(data, features = "AC013287.1", order = T) # location in the clusters

#Cluster3
FeaturePlot(data, features = c("ITM2C", "PIGR", "LINC01480", "MRM1", "LILRA4", "LILRB4", "MS4A1", "ISOC2", "MIR3681HG", "ALDOC"), order = T) # location in the clusters
FeaturePlot(data, features = c("ITM2C", "PIGR", "LINC01480", "LILRA4"), order = T) # location in the clusters

FeaturePlot(data, features = "ITM2C", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "PIGR", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "LINC01480", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "LILRA4", order = T, min.cutoff = 1) # location in the clusters

#Cluster4
FeaturePlot(data, features = c("LINC00987", "KLRC3", "AKR1C3", "IGFBP3", "HBA1", "AC007278.2", "TRGC1", "FCGR3A", "TRGV9", "AC130469.1"), order = T) # location in the clusters
FeaturePlot(data, features = c("FCGR3A", "KLRC3", "AKR1C3", "HBA1"), order = T) # location in the clusters

FeaturePlot(data, features = "FCGR3A", order = T) # location in the clusters
FeaturePlot(data, features = "KLRC3", order = T) # location in the clusters
FeaturePlot(data, features = "AKR1C3", order = T) # location in the clusters
FeaturePlot(data, features = "HBA1", order = T) # location in the clusters

#Cluster5
FeaturePlot(data, features = c("ADAM23", "RBMS3-AS2", "AC009041.2", "SPRY1", "MAL", "AC105402.3", "FAM153B", "NEFL", "MB21D2", "KRT72" ), order = T) # location in the clusters
FeaturePlot(data, features = c("ADAM23", "RBMS3-AS2", "MAL", "FAM153B"), order = T) # location in the clusters

FeaturePlot(data, features = "ADAM23", order = T) # location in the clusters
FeaturePlot(data, features = "RBMS3-AS2", order = T) # location in the clusters
FeaturePlot(data, features = "MAL", order = T) # location in the clusters
FeaturePlot(data, features = "FAM153B", order = T) # location in the clusters

#Cluster6

FeaturePlot(data, features = c("FCN1", "PADI4", "PID1", "ASGR2", "MAFB", "CCR1", "FAM20C", "LINC02202", "AC091138.1", "TLR5" ), order = T) # location in the clusters
FeaturePlot(data, features = c("FCN1"), order = T) # location in the clusters

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

grep('HBA', mark$cluster0)
grep('HBB', mark$cluster0)

# Subsclustering 

Idents(S2) <- S2$seurat_clusters
S2 <- FindNeighbors(S2, dims = 1:20, graph.name = "graph") #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.
S2 <- FindSubCluster(S2, resolution = 0.5, cluster = '4', graph.name = "graph")
table(S2$sub.cluster)

Idents(S2) <- as.factor(S2$sub.cluster)

DimPlot(S2, reduction = "umap", label = TRUE, pt.size = 0.5,  cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','#440154ff', 'lightblue', 'violet', 'blue3', 'orange')) + NoLegend()

markers2 <- FindAllMarkers(S2, only.pos = T)
top2 <- cut_markers(clusters = levels(markers2$cluster),markers2, ntop=50)

write.csv(top2, 'top_markers_subclustering_S2.csv')

### CELLTYPE ANNOTATION (AFTER SEARCHING FOR SPECIFIC MARKERS IN CELLTYPIST)

S2 <- SetIdent(S2, value = as.factor(S2$sub.cluster))

cluster_annotation <- c('Leukemic cells 1', 'Leukemic cells 2', 'Leukemic cells 3', 'Leukemic cells 4 (MS4A1+)', 'CD8+ cytotoxic T cells', 'T cells (FCGR3A+)', 'T cells (GNLY+)', 'CD4+ T cells', 'Leukemic cells 5', 'Myeloid cells')
cell_type <- cluster_annotation[as.factor(levels(S2))]
S2$cell_type <- cell_type[as.factor(S2$sub.cluster)]

names(cluster_annotation) <- levels(S2)

S2 <- SetIdent(S2, value = as.factor(S2$cell_type))

DimPlot(S2, reduction = "umap", label = F, pt.size = 0.5, cols = c('violet','pink','mediumseagreen','mediumaquamarine','darkseagreen','seagreen', 'green4', '#F4D35E','palevioletred3', 'palevioletred1')) + NoLegend()
DimPlot(S2, reduction = "umap", label = F, pt.size = 0.5, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','#440154ff', 'lightblue', 'violet', 'blue3', 'orange')) + NoLegend()

## SHINY APP

library(ShinyCell)

setwd("/home/adeluca/Documents/")
scConf = createConfig(S2)
makeShinyApp(S2, scConf, gene.mapping = TRUE,
             shiny.title = "Sample2")
