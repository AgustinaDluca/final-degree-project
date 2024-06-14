---
title: "patient S_180574"
author: "Agustina De Luca"
date: "2024-01-15"
output: html_document
---
## SAMPLE1

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
library(grDevices)

# What we did was to first pass all the raw data into the cell ranger ARC in order to get the data we can work with. Once this has finished in the cluster I 
# downloaded the folder named 'filtered_feature_bc_matrix' where all the three different matrix outputs are located
setwd("/Users/adeluca/Documents/Documents/")

# Load the Seurat object back into R
S1 <- readRDS("Documents/Documents/scRNA_samples/S_180574/original_object/Seuratobject_S180574.rds")

### SAVE ###
# Assuming your Seurat object is named 'your_seurat_object'
saveRDS(S1, file = "S1.rds")

### How to load multiome data

counts <- Read10X_h5("Documents/Documents/scRNA_samples/S_180574/filtered_feature_bc_matrix.h5")
names(counts)

data <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# Bash
# tabix -p bed /home/adeluca/Documents/S_180574/atac_fragments.tsv.gz

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#seqlevelsStyle(annotation) <- "UCSC"
ucsc.levels <- stringr::str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels
genome(annotation) <- "hg38"

data[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = "Documents/Documents/scRNA_samples/S_180574/atac_fragments.tsv.gz",
  annotation = annotation
)

## I checked that the the ZAP70 was present in the 'annotation' genome we selected (hg38)
## https://www.ncbi.nlm.nih.gov/datasets/gene/id/dc117a1eacd2c3839fde463efbb051fc/

dim(data)
# 36601  7125
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

hist(S1$percent.mt, main = 'Distribution of mitochondrial genes in S1')
hist(S1$nFeature_RNA, main = 'Distribution of the nFeatures in S1')

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

### SAVE ###
# Visualize the different metrics with violin plots
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

# Visualize feature-feature relationship with scatterplots
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

summary(data@meta.data$nFeature_RNA) # It is really interesting to look at the summary of the number of features as it helps us to choose a threshold to remove problematic cells
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 501     921    1068    1253    1287    8138 

summary(data@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   2.703   3.528   3.796   4.543   9.994 

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
#  0    1    2    3    4    5 
# 2170 2047  805  468  423   54 

'''
Con nfeatures = 800
   0    1    2    3    4    5 
1954 1763 1111  657  428   54 
'''

head(Idents(data), 5) #Cluster assignments for each cell after the clustering process.

## 8. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

data <- RunUMAP(data, dims = 1:20)
### SAVE ###

png('/home/adeluca/Documents/S_180574/nCount_RNA.png')
DimPlot(data, reduction = 'umap', label = T, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','#440154ff'))
dev.off()

FeaturePlot(data, features = 'nCount_RNA')
FeaturePlot(data, features = 'nFeature_RNA')
FeaturePlot(data, features = 'percent.mt')
dim(data)


### SAVE ###
VlnPlot(data,features = "nFeature_RNA",log = T, sort = T)

## 9. FINDING DIFFERENTIALLY EXPRESSED FEATURES (cluster biomarkers)

markers <- FindAllMarkers(data, only.pos = T)
### SAVE ###
# Save my markers matrix
#write.csv(markers, file = "/home/adeluca/Documents/S_180574/markers_S180574.csv", row.names = T)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

DoHeatmap(data, features = top10$gene) + NoLegend()

library(matchSCore2)
top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=50)
top

### SAVE ###
#write.csv(top, file = "/home/adeluca/Documents/S_180574/top50markers_S180574.csv", row.names = T)

### SAVE ###
#Most common markers
FeaturePlot(data, features = c("MS4A1","GNLY", "CD3E", "CD14", "FCE51A", "FCGR3A", "LYZ", "PPBP", "CD8A"),order = T) # location in the clusters
VlnPlot(data, features = c("MS4A1","GNLY", "CD3E", "CD14", "FCE51A", "FCGR3A", "LYZ", "PPBP", "CD8A"), slot = "counts", log = TRUE) # RNA counts
### SAVE ###

mark <- read.csv('/home/adeluca/Documents/S_180574/markers/top_markers_S180574.csv')

#Cluster0
FeaturePlot(data, features = c("OTOGL", "AC002377.1", "PTPRG", "ABCA9-AS1", "ADTRP", "SRGAP3", "HRK", "SRP14-AS1", "SLC14A1", "ABCA9"), order = T) # location in the clusters

FeaturePlot(data, features = c("AC002377.1", "ADTRP", "PTPRG", "HRK"), order = T) # location in the clusters

FeaturePlot(data, features = "AC002377.1", order = T) # location in the clusters
FeaturePlot(data, features = "ADTRP", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "HRK", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "PTPRG", order = T, min.cutoff = 1) # location in the clusters

#Cluster1
png('/home/adeluca/Documents/S_180574/clusters/cluster1_S1.png')
FeaturePlot(data, features = c("MAP2", "NTNG1", "SLC4A8", "IL13RA1", "LINC02427", "ARHGAP22", "IFNLR1", "FCRL5", "ARHGAP6", "LINC01480"), order = T) # location in the clusters
FeaturePlot(data, features = "MAP2", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "NTNG1", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "IFNLR1", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "SLC4A8", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "ARHGAP6", order = T, min.cutoff = 1) # location in the clusters

FeaturePlot(data, features = c("MAP2", "NTNG1", "SLC4A8"), order = T) # location in the clusters

#Cluster2
FeaturePlot(data, features = c("AL358115.1", "HMMR", "EXO1", "HGH1", "AL683807.2", "PCDHGA5", "RAB13", "AC087623.2", "AC092757.3", "PPM1E"), order = T) # location in the clusters
FeaturePlot(data, features = "AL358115.1", order = T) # location in the clusters
FeaturePlot(data, features = "HMMR", order = T) # location in the clusters
FeaturePlot(data, features = "EXO1", order = T) # location in the clusters
FeaturePlot(data, features = "HGH1", order = T) # location in the clusters

FeaturePlot(data, features = c("AL358115.1", "HMMR", "EXO1", "HGH1"), order = T) # location in the clusters

#Cluster3
FeaturePlot(data, features = c("AC005291.2", "RPS4Y2", "GARS", "RNF135", "UBE2Z", "BCAP29", "MTHFD2L", "ZSCAN16-AS1", "SYMPK", "ITFG2"), order = T) # location in the clusters
FeaturePlot(data, features = "AC005291.2", order = T, min.cutoff = 1) # location in the clusters
FeaturePlot(data, features = "RPS4Y2", order = T, min.cutoff = 1) # location in the clusters

FeaturePlot(data, features = c("AC005291.2", "RPS4Y2", "GARS"), order = T) # location in the clusters

#Cluster4
FeaturePlot(data, features = c("STAT4", "AC243829.2", "NCAM1", "AC010609.1", "PTGDS", "SH2D1A", "TNIP3", "SH2D2A", "AREG", "CD8A"), order = T) # location in the clusters
FeaturePlot(data, features = "STAT4", order = T) # location in the clusters
FeaturePlot(data, features = "SH2D1A", order = T) # location in the clusters
FeaturePlot(data, features = "AREG", order = T) # location in the clusters
FeaturePlot(data, features = "CD8A", order = T) # location in the clusters

FeaturePlot(data, features = c("STAT4", "SH2D1A", "AREG", "CD8A"), order = T) # location in the clusters

#Cluster5
FeaturePlot(data, features = c("TNFAIP2", "CD14", "CSF3R", "CLEC7A", "TIMP2", "TMTC1", "STAB1", "LRRC25", "SH3RF1", "NLRP3" ), order = T) # location in the clusters
FeaturePlot(data, features = "TNFAIP2", order = T) # location in the clusters
FeaturePlot(data, features = "CSF3R", order = T) # location in the clusters
FeaturePlot(data, features = "CLEC7A", order = T) # location in the clusters
FeaturePlot(data, features = "CD14", order = T) # location in the clusters

FeaturePlot(data, features = c("TNFAIP2", "CD14", "CSF3R", "CLEC7A"), order = T, min.cutoff = 1) # location in the clusters

# Color each cluster
DimPlot(data, group.by = "seurat_clusters", cols = c('#118AB2','grey', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','#06D6A0', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', '#8338EC', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', '#F4D35E', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', '#FF6B6B', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'grey', '#440154ff'))

grep('IL7R', mark$cluster0) # 4
grep('CCR7', mark$cluster0)
grep('S100A4', mark$cluster0) #2
grep('CD14', mark$cluster0)
grep('LYZ', mark$cluster0)
grep('MS4A1', mark$cluster0)
grep('CD8A', mark$cluster0)
grep('FCGR3A', mark$cluster0)
grep('MS4A7', mark$cluster0)
grep('GNLY', mark$cluster0) # 4
grep('NKG7', mark$cluster0) # 4
grep('PPBP', mark$cluster0)
grep('FCER1A', mark$cluster0)
grep('CST3', mark$cluster0)


# Subsclustering 
Idents(S1) <- S1$seurat_clusters
S1 <- FindNeighbors(S1, dims = 1:20, graph.name = "graph") #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.
S1 <- FindSubCluster(S1, resolution = 0.2, cluster = '4', graph.name = "graph")
table(S1$sub.cluster)

Idents(S1) <- as.factor(S1$sub.cluster)

DimPlot(S1, reduction = "umap", label = TRUE, pt.size = 0.5,  cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','#440154ff', 'red', 'violet')) + NoLegend()

markers1 <- FindAllMarkers(S1, only.pos = T)
top1 <- cut_markers(clusters = levels(markers1$cluster),markers1, ntop=50)

write.csv(top1, 'top_markers_subclustering_S1.csv')

### CELLTYPE ANNOTATION (AFTER SEARCHING FOR SPECIFIC MARKERS IN CELLTYPIST)

S1 <- SetIdent(S1, value = as.factor(S1$sub.cluster))
cluster_annotation <- c("Leukemic cells 1", "Leukemic cells 2 (MS4A1+)", "Leukemic cells 3", "Leukemic cells 4", "CD8+ Cytotoxic T cells","CD4+ T cells", "Cytotoxic NK", "CD14+ Mono")

cell_type <- cluster_annotation[as.factor(levels(S1))]
S1$cell_type <- cell_type[as.factor(S1$sub.cluster)]

names(cluster_annotation) <- levels(S1)

S1 <- SetIdent(S1, value = as.factor(S1$cell_type))

DimPlot(S1, reduction = "umap", label = F, pt.size = 0.5,  cols = c('#F4D35E','violet','pink','palevioletred','mediumseagreen','mediumaquamarine','darkseagreen','seagreen')) + NoLegend()

DimPlot(S1, reduction = "umap", label = F, pt.size = 0.5,  ccols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','lightblue', 'violet', '#440154ff')) + NoLegend()

## SHINY APP
library(ShinyCell)

setwd("/home/adeluca/Documents/")
scConf = createConfig(S1)
makeShinyApp(S1, scConf, gene.mapping = TRUE,
             shiny.title = "Sample1")
