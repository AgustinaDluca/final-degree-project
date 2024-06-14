---
title: "patient S_211855"
author: "Agustina De Luca"
date: "2024-01-15"
output: html_document
---

## SAMPLE4

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
S4 <- readRDS("S4.rds")

# Assuming your Seurat object is named 'your_seurat_object'
saveRDS(S4, file = "S4.rds")

### SAVE ###
# Assuming your Seurat object is named 'your_seurat_object'
#saveRDS(data, file = "/home/adeluca/Documents/S_211855/data_S211855.rds")
#saveRDS(data, file = "/home/adeluca/Documents/S_211855/Seuratobject_S211855.rds")

### How to load multiome data

counts <- Read10X_h5("/home/adeluca/Documents/S_211855/filtered_feature_bc_matrix.h5")
names(counts)

data <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# Bash
# tabix -p bed /home/adeluca/Documents/S_211855/atac_fragments.tsv.gz

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#seqlevelsStyle(annotation) <- "UCSC"
ucsc.levels <- stringr::str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels
genome(annotation) <- "hg38"

data[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = "/home/adeluca/Documents/S_211855/atac_fragments.tsv.gz",
  annotation = annotation
)

## I checked that the the ZAP70 was present in the 'annotation' genome we selected (hg38)
## https://www.ncbi.nlm.nih.gov/datasets/gene/id/dc117a1eacd2c3839fde463efbb051fc/

dim(data)
# 36601  7077
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
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

hist(S4$percent.mt, main = 'Distribution of mitochondrial genes in S4')
hist(S4$nFeature_RNA, main = 'Distribution of the nFeatures in S4')

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

data <- subset(data, subset = nFeature_RNA > 800 & percent.mt < 10) # We are going to select the cells we want to be analysed based on the previous results mainly the violin plots from above

# I could put a threshold in the percent.mt in 5 however it is not needed as the values as low it doesn't overcome the more than the 20 percent of the samples

### SAVE ###
# Visualize the different metrics with violin plots
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

# Visualize feature-feature relationship with scatterplots
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

summary(data@meta.data$nFeature_RNA) # It is really interesting to look at the summary of the number of features as it helps us to choose a threshold to remove problematic cells
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   809    1243    1412    1611    1684    8918 

summary(data@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2697  4.3312  5.5556  5.7388  7.0505  9.9877 

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

## 7. CLUSTER THE CELLS

data <- FindNeighbors(data, dims = 1:20) #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.

data <- FindClusters(data, resolution = 0.5) #Partition the KNN graph into clusters.

table(data@active.ident)
#    0    1    2    3    4 
#  3543 2002  367  287  132 

head(Idents(data), 5) #Cluster assignments for each cell after the clustering process.

## 8. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

data <- RunUMAP(data, dims = 1:20)
### SAVE ###
DimPlot(data, reduction = 'umap', label = T, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B'))

FeaturePlot(data, features = 'nCount_RNA')
FeaturePlot(data, features = 'nFeature_RNA')
FeaturePlot(data, features = 'percent.mt')

### SAVE ###
VlnPlot(data,features = "nFeature_RNA",log = T, sort = T)

## 9. FINDING DIFFERENTIALLY EXPRESSED FEATURES (cluster biomarkers)

markers <- FindAllMarkers(data, only.pos = T)
### SAVE ###
# Save my markers matrix
#write.csv(markers, file = "/home/adeluca/Documents/S_211855/markers_211855.csv", row.names = T)


markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(data, features = top10$gene) + NoLegend()

library(matchSCore2)
top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=50)
top
### SAVE ###
#write.csv(top, file = "/home/adeluca/Documents/S_211855/top50markers_S211855.csv", row.names = T)

### SAVE ###
#Most common markers
FeaturePlot(data, features = c("MS4A1","GNLY", "CD3E", "CD14", "FCE51A", "FCGR3A", "LYZ", "PPBP", "CD8A"),order = T) # location in the clusters
VlnPlot(data, features = c("MS4A1","GNLY", "CD3E", "CD14", "FCE51A", "FCGR3A", "LYZ", "PPBP", "CD8A"), slot = "counts", log = TRUE) # RNA counts

mark4 <- read.csv('S_211855/markers/top_markers_S211855.csv')

#Cluster0
FeaturePlot(data, features = c("CXCR4", "AL031846.2", "SAP25", "MED16", "UBR7", "KCTD5", "AC004158.1", "AC010978.1", "TIMM22", "SNAPIN"), order = T, min.cutoff = 1) # location in the clusters

FeaturePlot(data, features = c("SAP25", "MED16", "AC010978.1", "AC004158.1"), order = T) # location in the clusters
FeaturePlot(data, features = "SAP25", order = T) # location in the clusters
FeaturePlot(data, features = "MED16", order = T) # location in the clusters
FeaturePlot(data, features = "AC010978.1", order = T) # location in the clusters
FeaturePlot(data, features = "AC004158.1", order = T) # location in the clusters

#Cluster1
FeaturePlot(data, features = c("AL591845.1", "UNC79", "AC008691.1", "TNKS1BP1", "AC116407.2", "LINC01730", "LINC01833", "TAS2R30", "AC009065.4", "AC091057.2"), order = T) # location in the clusters

FeaturePlot(data, features = c("UNC79", "TNKS1BP1", "HIST1H2BB", "FHL1"), order = T) # location in the clusters
FeaturePlot(data, features = "UNC79", order = T) # location in the clusters
FeaturePlot(data, features = "TNKS1BP1", order = T) # location in the clusters
FeaturePlot(data, features = "LINC01730", order = T) # location in the clusters
FeaturePlot(data, features = "LINC01833", order = T) # location in the clusters

#Cluster2
FeaturePlot(data, features = c("PYCR1", "AC022784.1", "PALD1", "AC107993.1", "LRRC32", "AL136304.1", "LINC01136", "EGR3", "STUM", "Z94057.1"), order = T) # location in the clusters

FeaturePlot(data, features = c("PYCR1", "STUM", "AC022784.1"), order = T) # location in the clusters
FeaturePlot(data, features = "PYCR1", order = T) # location in the clusters
FeaturePlot(data, features = "STUM", order = T) # location in the clusters
FeaturePlot(data, features = "AC022784.1", order = T) # location in the clusters

#Cluster3
FeaturePlot(data, features = c("AC026803.1", "LINC00652", "MAPK12", "ARNTL2", "LURAP1L", "AC008543.5", "NECAB1", "AC090679.2", "GLIS2", "ADORA2B"), order = T) # location in the clusters

FeaturePlot(data, features = c("MAPK12", "ARNTL2", "GLIS2", "ADORA2B"), order = T) # location in the clusters
FeaturePlot(data, features = "MAPK12", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "ARNTL2", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "GLIS2", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "ADORA2B", order = T, min.cutoff = 1) # location in the clusters

#Cluster4
FeaturePlot(data, features = c("GPRIN3", "STAT4", "KLRK1", "CERS6", "TNFAIP2", "THEMIS", "CD28", "GIMAP7", "EIF4E3", "GATA3", "IL7R"), order = T) # location in the clusters

FeaturePlot(data, features = c("STAT4", "TNFAIP2", "THEMIS", "CD28", "IL7R"), order = T) # location in the clusters
FeaturePlot(data, features = "STAT4", order = T) # location in the clusters
FeaturePlot(data, features = "TNFAIP2", order = T) # location in the clusters
FeaturePlot(data, features = "THEMIS", order = T) # location in the clusters
FeaturePlot(data, features = "CD28", order = T) # location in the clusters
FeaturePlot(data, features = "IL7R", order = T) # location in the clusters


# Color each cluster

DimPlot(data, group.by = "seurat_clusters", cols = c('#118AB2','grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','#06D6A0', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', '#8338EC', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', '#F4D35E', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', '#FF6B6B'))

grep('IL7R', mark$cluster4)
grep('CCR7', mark$cluster4)
grep('S100A4', mark$cluster4)
grep('CD14', mark$cluster4)
grep('LYZ', mark$cluster4)
grep('MS4A1', mark$cluster4)
grep('CD8A', mark$cluster4)
grep('FCGR3A', mark$cluster4)
grep('MS4A7', mark$cluster4)
grep('GNLY', mark$cluster4)
grep('NKG7', mark$cluster4)
grep('PPBP', mark$cluster4)
grep('FCER1A', mark$cluster4)
grep('CST3', mark$cluster4)
grep('HBA', mark$cluster4)
grep('HBB', mark$cluster4)

# Subsclustering 
# Step 1: Perform subclustering
Idents(S4) <- S4$seurat_clusters
S4 <- FindNeighbors(S4, dims = 1:20, graph.name = "graph") #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.
S4 <- FindSubCluster(S4, cluster = '4', graph.name = "graph")

table(S4$sub.cluster)

Idents(S4) <- as.factor(S4$sub.cluster)

DimPlot(S4, reduction = "umap", label = TRUE, pt.size = 0.5,  cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','#440154ff', 'red', 'violet', 'blue', 'orange')) + NoLegend()

markers4 <- FindAllMarkers(S4, only.pos = T)
top4 <- cut_markers(clusters = levels(markers4$cluster),markers4, ntop=50)

write.csv(top4, 'top_markers_subclustering_S4.csv')

### CELLTYPE ANNOTATION (AFTER SEARCHING FOR SPECIFIC MARKERS IN CELLTYPIST)

S4 <- SetIdent(S4, value = as.factor(S4$sub.cluster))

cluster_annotation <- c("Leukemic cells 2", "Leukemic cells 1", "Leukemic cells 3", "CD4+ T cells", "CD8+ Cytotoxic T cells", "Leukemic cells 4")
cell_type <- cluster_annotation[as.factor(levels(S4))]
S4$cell_type <- cell_type[as.factor(S4$sub.cluster)]

names(cluster_annotation) <- levels(S4)

S4 <- SetIdent(S4, value = as.factor(S4$cell_type))

DimPlot(S4, reduction = "umap", label = TRUE, pt.size = 0.5, cols = c('#06D6A0','#118AB2','#8338EC','#FF6B6B', 'violet', '#F4D35E')) + NoLegend()
DimPlot(S4, reduction = "umap", label = F, pt.size = 0.5, cols = c('violet','pink','mediumseagreen','mediumaquamarine','darkseagreen','seagreen')) + NoLegend()

## SHINY APP
setwd("/home/adeluca/Documents/")
scConf = createConfig(S4)
makeShinyApp(S4, scConf, gene.mapping = TRUE,
             shiny.title = "Sample4")

# SIMULACION - Azimuth Seurat











