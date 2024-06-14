---
title: "all patients"
author: "Agustina De Luca"
date: "2024-01-23"
output: html_document
---
library(dplyr)
library(Seurat)
library(patchwork)
library(devtools)
library(hdf5r)
library(Signac)
library(biovizBase)
library(matchSCore2)
####################################

## MARKERS
markers_S180574 <- read.csv("../adeluca/Documents/S_180574/results_S180574/markers/markers_S180574.csv", header = F, row.names = 1)
markers_S210183 <- read.csv("../adeluca/Documents/S_210183/results_S210183/markers/markers_S210183.csv", header = F, row.names = 1)
markers_S210414 <- read.csv("../adeluca/Documents/S_210414/results_S210414/markers/markers_S210414.csv", header = F, row.names = 1)
markers_S211855 <- read.csv("../adeluca/Documents/S_211855/results_S211855/makers/markers_211855.csv", header = F, row.names = 1)

## TOP MARKERS
top_markers_S180574 <- read.csv("../adeluca/Documents/top_markers_S180574.csv")
top_markers_S210183 <- read.csv("../adeluca/Documents/top_markers_S210183.csv")
top_markers_S210414 <- read.csv("../adeluca/Documents/top_markers_S210414.csv")
top_markers_S211855 <- read.csv("../adeluca/Documents/top_markers_S211855.csv")

S1 <- readRDS("../adeluca/Documents/S1.rds")
# UMAP plot
DimPlot(S1, reduction = 'umap', label = T, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','lightblue', 'violet', '#440154ff'))

# Vln plot of all clusters
VlnPlot(S1,features = "nFeature_RNA",log = T, sort = T, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','lightblue', 'violet', '#440154ff'))

S2 <- readRDS("../adeluca/Documents/S2.rds")
# UMAP plot
DimPlot(S2, reduction = 'umap', label = T, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','lightblue', 'violet', '#440154ff'))
# Vln plot of all clusters
VlnPlot(S2,features = "nFeature_RNA",log = T, sort = T, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','#440154ff', 'violet', 'lightblue', 'darkgreen', 'deepskyblue2'))

S3 <- readRDS("../adeluca/Documents/S3.rds")
# UMAP plot
DimPlot(S3, reduction = 'umap', label = T, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','#440154ff', 'violet'))
# Vln plot of all clusters
VlnPlot(S3,features = "nFeature_RNA",log = T, sort = T, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','#440154ff', 'violet'))

S4 <- readRDS("S4.rds")
# UMAP plot
DimPlot(S4, reduction = 'umap', label = T, cols = c('#118AB2','#06D6A0','#8338EC','#FF6B6B', 'violet', '#F4D35E'))
# Vln plot of all clusters
VlnPlot(S4,features = "nFeature_RNA",log = T, sort = T, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B', '#440154ff'))

#Most common markers
FeaturePlot(data, features = c("MS4A1","GNLY", "CD3E", "CD14", "FCE51A", "FCGR3A", "LYZ", "PPBP", "CD8A"),order = T) # location in the clusters
VlnPlot(data, features = c("MS4A1","GNLY", "CD3E", "CD14", "FCE51A", "FCGR3A", "LYZ", "PPBP", "CD8A"), slot = "counts", log = TRUE) # RNA counts

head(top_markers_S180574, 10)
head(top_markers_S210183, 10)
head(top_markers_S210414, 10)
head(top_markers_S211855, 10)

#cluster5.markers <- FindMarkers(S_180574, ident.1 = 5, ident.2 = c(0, 3))
#head(cluster5.markers, n = 5)

FeaturePlot(S1, features = 'ZAP70',order = T) # location in the clusters
FeaturePlot(S2, features = 'ZAP70',order = T) # location in the clusters
FeaturePlot(S3, features = 'ZAP70',order = T) # location in the clusters
FeaturePlot(S4, features = 'ZAP70',order = T) # location in the clusters

VlnPlot(S1, features = 'ZAP70', slot = "counts", log = TRUE, cols = c('#118AB2', '#06D6A0', '#F4D35E', '#FF6B6B', '#8338EC','#440154ff','#FF6B6B','#8338EC', '#06D6A0')) # RNA counts
VlnPlot(S2, features = 'ZAP70', slot = "counts", log = TRUE, cols = c('#118AB2', '#F4D35E', '#06D6A0', '#8338EC', '#06D6A0','#FF6B6B','#FF6B6B','#440154ff', '#F4D35E', 'violet', '#8338EC', '#440154ff', '#118AB2', 'violet')) # RNA counts
VlnPlot(S3, features = 'ZAP70', slot = "counts", log = TRUE, cols = c('#8338EC', '#06D6A0', '#F4D35E', '#118AB2', '#06D6A0','#FF6B6B','#118AB2','#FF6B6B', '#F4D35E', '#440154ff', '#8338EC')) # RNA counts
VlnPlot(S4, features = 'ZAP70', slot = "counts", log = TRUE, cols = c('#118AB2', '#06D6A0', '#8338EC', '#F4D35E', '#FF6B6B')) # RNA counts

# Calculate the UMAP reduction if not already done
S1 <- RunUMAP(S1, reduction = "umap")

# Create the plot with the specified gene expression
DimPlot(S1, group.by = 'ZAP70', label = F)
DimPlot(S2, group.by = 'ZAP70', label = F)
DimPlot(S3, group.by = 'ZAP70', label = F)
DimPlot(S4, group.by = 'ZAP70', label = F)

## Label cells that express the ZAP70 marker in more than one count per cell

# Create a column in the metadata for ZAP expression, here we will evaluate if a cell is ZAP positive or not

S1@meta.data$ZAP_expression <- ifelse(S1@assays$RNA$counts['ZAP70', ] >= 1, "ZAP+", "ZAP-")

table(S1@meta.data$ZAP_expression)

#Load the object
S1 <- readRDS("../adeluca/Documents/S_180574/data_S180574.rds")
#Add a column regarding the ZAP70 expression
S1@meta.data$ZAP_expression <- ifelse(S1@assays$RNA$counts['ZAP70', ] >= 1, "ZAP+", "ZAP-")
#Look at the proportion
table(S1@meta.data$ZAP_expression)
# ZAP- ZAP+ 
# 5904   63 
saveRDS(S1, file = "/home/adeluca/Documents/S1.rds")


#Load the object
S2 <- readRDS("../adeluca/Documents/S_210183/data_S210183.rds")
#Add a column regarding the ZAP70 expression
S2@meta.data$ZAP_expression <- ifelse(S2@assays$RNA$counts['ZAP70', ] >= 1, "ZAP+", "ZAP-")
#Look at the proportion
table(S2@meta.data$ZAP_expression)
# ZAP- ZAP+ 
# 6912  908 
saveRDS(S2, file = "/home/adeluca/Documents/S2.rds")


#Load the object
S3 <- readRDS("../adeluca/Documents/S_210414/data_S210414.rds")
#Add a column regarding the ZAP70 expression
S3@meta.data$ZAP_expression <- ifelse(S3@assays$RNA$counts['ZAP70', ] >= 1, "ZAP+", "ZAP-")
#Look at the proportion
table(S3@meta.data$ZAP_expression)
# ZAP- ZAP+ 
# 838  236 
saveRDS(S3, file = "new_S3/S3_new.rds")


#Load the object
S4 <- readRDS("../adeluca/Documents/S_211855/data_S211855.rds")
#Add a column regarding the ZAP70 expression
S4@meta.data$ZAP_expression <- ifelse(S4@assays$RNA$counts['ZAP70', ] >= 1, "ZAP+", "ZAP-")
#Look at the proportion
table(S4@meta.data$ZAP_expression)
# ZAP- ZAP+ 
# 6236   95
saveRDS(S4, file = "/home/adeluca/Documents/S4.rds")

S1 <- readRDS("/home/adeluca/Documents/S1.rds")
S2 <- readRDS("/home/adeluca/Documents/S2.rds")
S3 <- readRDS("/home/adeluca/Documents/S3.rds")
S4 <- readRDS("/home/adeluca/Documents/S4.rds")

table(S1@meta.data$seurat_clusters, S1@meta.data$ZAP_expression)
#     ZAP- ZAP+
# 0  2170    0
# 1  2046    1
# 2  801     4
# 3  468     0
# 4  365    58
# 5   54     0
table(S2@meta.data$seurat_clusters, S2@meta.data$ZAP_expression)
#    ZAP- ZAP+
# 0  2143  121
# 1  1875  146
# 2  1168  283
# 3  1175  106
# 4  334  122
# 5  195  127
# 6   22    3
table(S3@meta.data$seurat_clusters, S3@meta.data$ZAP_expression)
#    ZAP- ZAP+
# 0  326   45
# 1  196   21
# 2  184   17
# 3   57  109
# 4   59   44
# 5   16    0
table(S4@meta.data$seurat_clusters, S4@meta.data$ZAP_expression)
#    ZAP- ZAP+
# 0  3511   32
# 1  1980   22
# 2  365    2
# 3  283    4
# 4   97   35

#head(data@meta.data[, 8:9], 10)

## PLOT PER CLUSTER ZAP+/ZAP-
library(ggplot2)

proportions <- as.data.frame((table(S4$seurat_clusters, S4$ZAP_expression) / rowSums(table(S4$seurat_clusters, S4$ZAP_expression)))*100)

# Assuming 'proportions' is your S1 frame
ZAP70_expression_plot <- ggplot(proportions, aes(x = factor(Var1), y = Freq, fill = Var2)) +
                          geom_bar(stat = "identity", position = "stack") +
                          scale_y_continuous(labels = scales::percent_format(scale = 1)) +
                          labs(title = "Proportion of ZAP+ and ZAP- in each cluster in S4",
                               x = "Clusters",
                               y = "Percentage",
                               fill = "Expression") +
                          theme_classic() +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
                          scale_fill_manual(values = c("ZAP-" = "skyblue", "ZAP+" = "salmon"))

ggsave("ZAP70_expressionS4.png", ZAP70_expression_plot, width = 6, height = 4)

## UMAP - ZAP70+

S1@meta.data$ZAP_expression <- as.factor(S1@meta.data$ZAP_expression)
S2@meta.data$ZAP_expression <- as.factor(S2@meta.data$ZAP_expression)
S3@meta.data$ZAP_expression <- as.factor(S3@meta.data$ZAP_expression)
S4@meta.data$ZAP_expression <- as.factor(S4@meta.data$ZAP_expression)

# EXPRESSION
plot <- FeaturePlot(S1, features = 'ZAP70', order = TRUE, cols = c("grey", "red"), pt.size = 0.5)
ggsave(filename = "/home/adeluca/Documents/umap_ZAP70+_S1.jpg", height = 7, width = 12, plot = plot, quality = 50)

plot <- FeaturePlot(S2, features = 'ZAP70', order = TRUE, cols = c("grey", "red"), pt.size = 0.5)
ggsave(filename = "/home/adeluca/Documents/umap_ZAP70+_S2.jpg", height = 7, width = 12, plot = plot, quality = 50)

plot <- FeaturePlot(S3, features = 'ZAP70', order = TRUE, cols = c("grey", "red"), pt.size = 0.5)
ggsave(filename = "/home/adeluca/Documents/umap_ZAP70+_S3.jpg", height = 7, width = 12, plot = plot, quality = 50)

plot <- FeaturePlot(S4, features = 'ZAP70', order = TRUE, cols = c("grey", "red"), pt.size = 0.5)
ggsave(filename = "/home/adeluca/Documents/umap_ZAP70+_S4.jpg", height = 7, width = 12, plot = plot, quality = 50)

# COLOR CELLS THAT ARE ZAP+ 
plot <- DimPlot(S1, group.by = "ZAP_expression", cols = c("lightgray", "red"), pt.size = 0.1, order = T)
ggsave(filename = "/home/adeluca/Documents/umap_ZAP70_S1.jpg", height = 7, width = 12, plot = plot, quality = 50)

plot <- DimPlot(S2, group.by = "ZAP_expression", cols = c("lightgray", "red"), pt.size = 0, order = T)
ggsave(filename = "/home/adeluca/Documents/umap_ZAP70_S2.jpg", height = 7, width = 12, plot = plot, quality = 50)

plot <- DimPlot(S3, group.by = "ZAP_expression", cols = c("lightgray", "red"), pt.size = 0.1, order = T)
ggsave(filename = "/home/adeluca/Documents/umap_ZAP70_S3.jpg", height = 7, width = 12, plot = plot, quality = 50)

plot <- DimPlot(S4, group.by = "ZAP_expression", cols = c("lightgray", "red"), pt.size = 0.1, order = T)
ggsave(filename = "/home/adeluca/Documents/umap_ZAP70_S4.jpg", height = 7, width = 12, plot = plot, quality = 50)

##Create a violin plot for ZAP70 in ZAP70+ cells
#S_[, S_$ZAP_expression == 'ZAP+'] -> FILTER cells that are ZAP70+

VlnPlot(S1[, S1$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'sub.cluster', log = TRUE, sort = T, cols = c('#FF6B6B','violet','lightblue','#8338EC', '#06D6A0'))#S1 - Groups with fewer than two data points have been dropped. 
VlnPlot(S2[, S2$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'sub.cluster', log = TRUE, sort = T, cols = c('#118AB2','#FF6B6B','lightblue','#06D6A0','#F4D35E','#8338EC', 'violet', '#440154ff'))
VlnPlot(S3[, S3$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'seurat_clusters', log = TRUE, sort = T, cols = c('#118AB2','#06D6A0','#8338EC','#F4D35E','#FF6B6B','#440154ff', 'violet'))
VlnPlot(S4[, S4$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'sub.cluster', log = TRUE, sort = T, cols = c('violet','#06D6A0','#8338EC','#FF6B6B', '#118AB2', '#F4D35E'))

# DE-GENES BETWEEN ZAP+ ANZ ZAP-

S1 <- SetIdent(S1, value = S1$cluster_expression)
DE_genes <- FindMarkers(object = S1, ident.1 = 'ZAP+', ident.2 = 'ZAP-') # Run differential expression analysis
#write.csv(DE_genes, file = "/home/adeluca/Documents/ZAP70/S1/DE_genes_S1.csv", row.names = T)

S2 <- SetIdent(S2, value = S2$ZAP_expression)
DE_genes <- FindMarkers(object = S2, ident.1 = 'ZAP+', ident.2 = 'ZAP-') # Run differential expression analysis
#write.csv(DE_genes, file = "/home/adeluca/Documents/ZAP70/S2/DE_genes_S2.csv", row.names = T)

S3 <- SetIdent(S3, value = S3$ZAP_expression)
DE_genes <- FindMarkers(object = S3, ident.1 = 'ZAP+', ident.2 = 'ZAP-') # Run differential expression analysis
#write.csv(DE_genes, file = "/home/adeluca/Documents/ZAP70/S3/DE_genes_S3.csv", row.names = T)

S4 <- SetIdent(S4, value = S4$)
DE_genes <- FindMarkers(object = S4, ident.1 = 'ZAP+', ident.2 = 'ZAP-', logfc.threshold = 0.1) # Run differential expression analysis
#write.csv(DE_genes, file = "/home/adeluca/Documents/ZAP70/S4/DE_genes_S4.csv", row.names = T)

### IMPORTANTEEE ###

##Dot plot of ZAP70+/- markers per cluster of sample 
# The size of the dot corresponds to the percentage of cells expressing the feature in each cluster. 
# The color represents the average expression level
#ZAP+/- clusters
S4$cluster_expression <- paste(S4$seurat_clusters, S4$ZAP_expression, sep = "_")

S3 <- SetIdent(S3, value = as.factor(S3$cluster_expression))

markers = read.csv('markers_pos_neg.csv')
order = markers$ZAP70pos_neg_markers

DoHeatmap(S2, features = order) + NoLegend()

#ZAP+
DotPlot(S1, features = c('KLRC4', 'NPDC1', 'MEX3D', 'KRT73', 'HOPX', 'JAKMIP1', 'B3GAT1', 'AC104078.1', 'LINC01980', 'LINC02611', 'AC097511.1')) + RotatedAxis()
DotPlot(S2, features = c('KLRC4', 'NPDC1', 'MEX3D', 'KRT73', 'HOPX', 'JAKMIP1', 'B3GAT1', 'AC104078.1', 'LINC01980', 'LINC02611', 'AC097511.1')) + RotatedAxis()
DotPlot(S3, features = c('KLRC4', 'NPDC1', 'MEX3D', 'KRT73', 'HOPX', 'JAKMIP1', 'B3GAT1', 'AC104078.1', 'LINC01980', 'LINC02611', 'AC097511.1')) + RotatedAxis()
DotPlot(S4, features = c('KLRC4', 'NPDC1', 'MEX3D', 'KRT73', 'HOPX', 'JAKMIP1', 'B3GAT1', 'AC104078.1', 'LINC01980', 'LINC02611', 'AC097511.1')) + RotatedAxis()
#ZAP-
DotPlot(S1, features = c('PLBD1', 'VCAN', 'MNDA', 'CHIT1', 'SRGAP3', 'BTBD8', 'ZBTB20-AS1')) + RotatedAxis()
DotPlot(S2, features = c('PLBD1', 'VCAN', 'MNDA', 'CHIT1', 'SRGAP3', 'BTBD8', 'ZBTB20-AS1')) + RotatedAxis()
DotPlot(S3, features = c('PLBD1', 'VCAN', 'MNDA', 'CHIT1', 'SRGAP3', 'BTBD8', 'ZBTB20-AS1')) + RotatedAxis()
DotPlot(S4, features = c('PLBD1', 'VCAN', 'MNDA', 'CHIT1', 'SRGAP3', 'BTBD8', 'ZBTB20-AS1')) + RotatedAxis()

##Dot plot of markers per cluster of each sample (agarre los primeros dos markers de cada clsuter)
DotPlot(S1, features = c('PTPRG', 'OTOGL', 'SLC4A8', 'NTNG1', 'EXO1', 'HGH1', 'AC005291.2', 'RPS4Y2', 'STAT4', 'CD8A', 'TNFAIP2', 'CD14'), col.min = 0.8, dot.scale = 5, scale.by = 'size') + RotatedAxis()
DotPlot(S2, features = c('LINC01239', 'AP003352.1', 'STK32B', 'SGSM1', 'POF1B', 'WNT6', 'ITM2C', 'PIGR', 'FCGR3A', 'KLRC3', 'ADAM23', 'MAL', 'FCN1', 'PADI4'), col.min = 0.1, dot.scale = 5, scale.by = 'size') + RotatedAxis()

DotPlot(S3, features = c("GAB1", "PAX5", "ACO2", "PIKFYVE", "CD3G", "STAT4", "CD19", "EGR1", "PLD5", "WNT2B", "VCAN", "TNFAIP2", "CADM2", "PTPRD"), col.min = 0.5, dot.scale = 5, scale.by = 'size') + RotatedAxis() # location in the clusters


DotPlot(S4, features = c('UBR7', 'KCTD5', 'AL591845.1', 'UNC79', 'PYCR1', 'STUM', 'GLIS2', 'ARNTL2', 'THEMIS', 'STAT4'), col.min = 0, dot.scale = 5, scale.by = 'size') + RotatedAxis()

##Dot plot of clusters/ZAP70 with common markers

DotPlot(S1, features = c("IL7R", "CCR7", "S100A4", "CD4","CD8A", "CD3E", "CD3G", "GNLY", "NKG7", "MS4A1", "FCGR3A", "MS4A7", "CD14", "LYZ", "FCER1A", "CST3", "PPBP"), col.min = 0) + RotatedAxis()
DotPlot(S2, features = c("IL7R", "CCR7", "S100A4", "CD4","CD8A", "CD3E", "CD3G", "GNLY", "NKG7", "MS4A1", "FCGR3A", "MS4A7", "CD14", "LYZ", "FCER1A", "CST3", "PPBP"), col.min = 0) + RotatedAxis()
DotPlot(S3, features = c("IL7R", "CCR7", "S100A4", "CD4","CD8A", "CD3E", "CD3G", "GNLY", "NKG7", "MS4A1", "FCGR3A", "MS4A7", "CD14", "LYZ", "FCER1A", "CST3", "PPBP"), col.min = 0) + RotatedAxis()
DotPlot(S4, features = c("IL7R", "CCR7", "S100A4", "CD4","CD8A", "CD3E", "CD3G", "GNLY", "NKG7", "MS4A1", "FCGR3A", "MS4A7", "CD14", "LYZ", "FCER1A", "CST3", "PPBP"), col.min = 0) + RotatedAxis()

## Heatmap w/ DE-genes from ZAP+ and ZAP-
## Get top 100 from ZAP+/ZAP-
# SAMPLE1
S1_DE_genes_ZAP70 <- read.csv('ZAP70/S1/DE_genes_S1.csv')

S1_ZAP70_pos <- head(S1_DE_genes_ZAP70$X, 500)
S1_ZAP70_neg <- tail(S1_DE_genes_ZAP70$X, 500)

S1_ZAP70 <- data.frame(S1_ZAP70_pos, S1_ZAP70_neg)
names(S1_ZAP70) <- c("ZAP+", "ZAP-")

# SAMPLE2
S2_DE_genes_ZAP70 <- read.csv('ZAP70/S2/DE_genes_S2.csv')

S2_ZAP70_pos <- head(S2_DE_genes_ZAP70$X, 500)
S2_ZAP70_neg <- tail(S2_DE_genes_ZAP70$X, 500)

S2_ZAP70 <- data.frame(S2_ZAP70_pos, S2_ZAP70_neg)
names(S2_ZAP70) <- c("ZAP+", "ZAP-")

# SAMPLE3
S3_DE_genes_ZAP70 <- read.csv('ZAP70/S3/DE_genes_S3.csv')

S3_ZAP70_pos <- head(S3_DE_genes_ZAP70$X, 500)
S3_ZAP70_neg <- tail(S3_DE_genes_ZAP70$X, 500)

S3_ZAP70 <- data.frame(S3_ZAP70_pos, S3_ZAP70_neg)
names(S3_ZAP70) <- c("ZAP+", "ZAP-")

# SAMPLE4
S4_DE_genes_ZAP70 <- read.csv('ZAP70/S4/DE_genes_S4.csv')

S4_ZAP70_pos <- head(S4_DE_genes_ZAP70$X, 500)
S4_ZAP70_neg <- tail(S4_DE_genes_ZAP70$X, 500)

S4_ZAP70 <- data.frame(S4_ZAP70_pos, S4_ZAP70_neg)
names(S4_ZAP70) <- c("ZAP+", "ZAP-")

library(matchSCore2)
## The matchSCore2 function computes the clustering comparison and produce the heatmap table with Jaccard Indexes for each group combination

ZAP70_12 <- matchSCore2(gene_cl.ref = S1_ZAP70, gene_cl.obs = S2_ZAP70, ylab = "SAMPLE1",xlab = "SAMPLE2")
ZAP70_13 <- matchSCore2(gene_cl.ref = S1_ZAP70, gene_cl.obs = S3_ZAP70, ylab = "SAMPLE1",xlab = "SAMPLE3")
ZAP70_14 <- matchSCore2(gene_cl.ref = S1_ZAP70, gene_cl.obs = S4_ZAP70, ylab = "SAMPLE1",xlab = "SAMPLE4")
ZAP70_23 <- matchSCore2(gene_cl.ref = S2_ZAP70, gene_cl.obs = S3_ZAP70, ylab = "SAMPLE2",xlab = "SAMPLE3")
ZAP70_24 <- matchSCore2(gene_cl.ref = S2_ZAP70, gene_cl.obs = S4_ZAP70, ylab = "SAMPLE2",xlab = "SAMPLE4")
ZAP70_34 <- matchSCore2(gene_cl.ref = S3_ZAP70, gene_cl.obs = S4_ZAP70, ylab = "SAMPLE3",xlab = "SAMPLE4")

## The matchSCore heatmap 
ZAP70_12$ggplot
ZAP70_13$ggplot
ZAP70_14$ggplot
ZAP70_23$ggplot
ZAP70_24$ggplot
ZAP70_34$ggplot

# ms$matchScore - average similarity (Jaccard Index), a higher matching score indicates a higher degree of similarity 
# ms$labels - labels associated with the maximum Jaccard Index values
# ms$max_JI - maximum Jaccard Index values for each reference group (TOP_1)
# ms$JI.mat - matrix represents the Jaccard Index values between all pairs of reference groups and cluster-specific gene sets 
# ms$ggplot - graphical representation of the comparison

## LOOK FOR SHARED MARKERS
library(VennDiagram)
library(scales)
library(gplots)

read.csv('S1_ZAP70_MARKERS.csv')
read.csv('S2_ZAP70_MARKERS.csv')
read.csv('S3_ZAP70_MARKERS.csv')
read.csv('S4_ZAP70_MARKERS.csv')

#replace each S1_ZAP70 por el nombre de la table

# ZAP70+

x = list(S1_ZAP70$`ZAP+`, S2_ZAP70$`ZAP+`, S3_ZAP70$`ZAP+`, S4_ZAP70$`ZAP+`)

venn.plot.ZAP70pos <- venn.diagram(
 x = list(S1_ZAP70$`ZAP+`, S2_ZAP70$`ZAP+`, S3_ZAP70$`ZAP+`, S4_ZAP70$`ZAP+`),
 category.names = c("Sample 1", "Sample 2", "Sample 3", "Sample 4"),
 filename = NULL,
 col = c('#FF6B6B', '#DAA520', '#440154ff', '#21908dff'),
 fill = c(alpha("#FF6B6B",0.3), alpha('#DAA520',0.3), alpha('#440154ff',0.3), alpha('#21908dff',0.3))
)
grid.draw(venn.plot.ZAP70pos)

sample1_2_pos <- intersect(x[[1]], x[[2]])
sample1_3_pos <- intersect(x[[1]], x[[3]])
sample1_4_pos <- intersect(x[[1]], x[[4]])
sample2_3_pos <- intersect(x[[2]], x[[3]])
sample2_4_pos <- intersect(x[[2]], x[[4]])
sample3_4_pos <- intersect(x[[3]], x[[4]])

sample1_2_3_pos <- Reduce(intersect, x[c(1, 2, 3)])
sample1_2_4_pos <- Reduce(intersect, x[c(1, 2, 4)])
sample2_3_4_pos <- Reduce(intersect, x[c(3, 2, 4)])
sample1_3_4_pos <- Reduce(intersect, x[c(1, 3, 4)])

all_samples_pos <- Reduce(intersect, x)

# Obtener la longitud máxima de todas las listas
max_length <- max(length(sample1_2_pos), length(sample1_3_pos), length(sample1_4_pos),
                  length(sample2_3_pos), length(sample2_4_pos), length(sample3_4_pos),
                  length(sample1_2_3_pos), length(sample1_2_4_pos), length(sample2_3_4_pos),
                  length(sample1_3_4_pos), length(all_samples_pos))

# Rellenar las listas más cortas con NA hasta que todas tengan la misma longitud
sample1_2_pos <- c(sample1_2_pos, rep(NA, max_length - length(sample1_2_pos)))
sample1_3_pos <- c(sample1_3_pos, rep(NA, max_length - length(sample1_3_pos)))
sample1_4_pos <- c(sample1_4_pos, rep(NA, max_length - length(sample1_4_pos)))
sample2_3_pos <- c(sample2_3_pos, rep(NA, max_length - length(sample2_3_pos)))
sample2_4_pos <- c(sample2_4_pos, rep(NA, max_length - length(sample2_4_pos)))
sample3_4_pos <- c(sample3_4_pos, rep(NA, max_length - length(sample3_4_pos)))
sample1_2_3_pos <- c(sample1_2_3_pos, rep(NA, max_length - length(sample1_2_3_pos)))
sample1_2_4_pos <- c(sample1_2_4_pos, rep(NA, max_length - length(sample1_2_4_pos)))
sample2_3_4_pos <- c(sample2_3_4_pos, rep(NA, max_length - length(sample2_3_4_pos)))
sample1_3_4_pos <- c(sample1_3_4_pos, rep(NA, max_length - length(sample1_3_4_pos)))
all_samples_pos <- c(all_samples_pos, rep(NA, max_length - length(all_samples_pos)))

# Fusionar todas las listas
sample_intersections <- Map(c, sample1_2_pos, sample1_3_pos, sample1_4_pos,
                            sample2_3_pos, sample2_4_pos, sample3_4_pos,
                            sample1_2_3_pos, sample1_2_4_pos, sample2_3_4_pos,
                            sample1_3_4_pos, all_samples_pos)

# Convertir la lista en un data frame
sample_intersections_df <- as.data.frame(sample_intersections, header = FALSE)
row.names(sample_intersections_df) <- c('sample1_2_pos', 'sample1_3_pos', 'sample1_4_pos', 'sample2_3_pos', 'sample2_4_pos', 'sample3_4_pos', 'sample1_2_3_pos', 'sample1_2_4_pos', 'sample2_3_4_pos', 'sample1_3_4_pos', 'all_samples_pos')

df <- t(sample_intersections_df)
rownames(df) <- NULL

write.csv(df, file = '/home/adeluca/Documents/venn_intersection_ZAP70pos.csv')

# ZAP70-

x = list(S1_ZAP70$`ZAP-`, S2_ZAP70$`ZAP-`, S3_ZAP70$`ZAP-`, S4_ZAP70$`ZAP-`)

venn.plot.ZAP70neg <- venn.diagram(
  x = list(S1_ZAP70$`ZAP-`, S2_ZAP70$`ZAP-`, S3_ZAP70$`ZAP-`, S4_ZAP70$`ZAP-`),
  category.names = c("Sample 1", "Sample 2", "Sample 3", "Sample 4"),
  filename = NULL,
  col = c('#FF6B6B', '#DAA520', '#440154ff', '#21908dff'),
  fill = c(alpha("#FF6B6B",0.3), alpha('#DAA520',0.3), alpha('#440154ff',0.3), alpha('#21908dff',0.3))
)
grid.draw(venn.plot.ZAP70neg)


sample1_2_neg <- intersect(x[[1]], x[[2]])
sample1_3_neg <- intersect(x[[1]], x[[3]])
sample1_4_neg <- intersect(x[[1]], x[[4]])
sample2_3_neg <- intersect(x[[2]], x[[3]])
sample2_4_neg <- intersect(x[[2]], x[[4]])
sample3_4_neg <- intersect(x[[3]], x[[4]])

sample2_3_4_neg <- Reduce(intersect, x[c(3, 2, 4)])
sample1_3_4_neg <- Reduce(intersect, x[c(1, 3, 4)])


# Obtener la longitud máxima de todas las listas
max_length <- max(length(sample1_2_neg), length(sample1_3_neg), length(sample1_4_neg),
                  length(sample2_3_neg), length(sample2_4_neg), length(sample3_4_neg),
                  length(sample2_3_4_neg), length(sample1_3_4_neg))

# Rellenar las listas más cortas con NA hasta que todas tengan la misma longitud
sample1_2_neg <- c(sample1_2_neg, rep(NA, max_length - length(sample1_2_neg)))
sample1_3_neg <- c(sample1_3_neg, rep(NA, max_length - length(sample1_3_neg)))
sample1_4_neg <- c(sample1_4_neg, rep(NA, max_length - length(sample1_4_neg)))
sample2_3_neg <- c(sample2_3_neg, rep(NA, max_length - length(sample2_3_neg)))
sample2_4_neg <- c(sample2_4_neg, rep(NA, max_length - length(sample2_4_neg)))
sample3_4_neg <- c(sample3_4_neg, rep(NA, max_length - length(sample3_4_neg)))
sample2_3_4_neg <- c(sample2_3_4_neg, rep(NA, max_length - length(sample2_3_4_neg)))
sample1_3_4_neg <- c(sample1_3_4_neg, rep(NA, max_length - length(sample1_3_4_neg)))

# Fusionar todas las listas
sample_intersections <- Map(c, sample1_2_neg, sample1_3_neg, sample1_4_neg,
                            sample2_3_neg, sample2_4_neg, sample3_4_neg,
                            sample2_3_4_neg, sample1_3_4_neg)

# Convertir la lista en un data frame
sample_intersections_df <- as.data.frame(sample_intersections, header = FALSE)
row.names(sample_intersections_df) <- c('sample1_2_neg', 'sample1_3_neg', 'sample1_4_neg', 'sample2_3_neg', 'sample2_4_neg', 'sample3_4_neg', 'sample2_3_4_neg', 'sample1_3_4_neg')

df <- t(sample_intersections_df)
rownames(df) <- NULL

write.csv(df, file = '/home/adeluca/Documents/venn_intersection_ZAP70neg.csv')

write.csv(markers, file = 'markers_pos_neg.csv')

markers  = as.data.frame(c(sample1_2_pos, sample1_3_pos, sample1_4_pos, sample2_3_pos, sample2_4_pos, sample3_4_pos, sample1_2_3_pos, sample1_2_4_pos, sample2_3_4_pos, sample1_3_4_pos, all_samples_pos, sample1_2_neg, sample1_3_neg, sample1_4_neg, sample2_3_neg, sample2_4_neg, sample3_4_neg, sample2_3_4_neg, sample1_3_4_neg))
names(markers) = 'ZAP70pos_neg_markers'

'''
pos_table <- venn(list(S1_ZAP70$`ZAP+`, S2_ZAP70$`ZAP+`, S3_ZAP70$`ZAP+`, S4_ZAP70$`ZAP+`))
neg_table <- venn(list(S1_ZAP70$`ZAP-`, S2_ZAP70$`ZAP-`, S3_ZAP70$`ZAP-`, S4_ZAP70$`ZAP-`))

saveRDS(pos_table, file = "/home/adeluca/Documents/ZAP70/ZAP+_table_all.samples.rds")
saveRDS(neg_table, file = "/home/adeluca/Documents/ZAP70/ZAP-_table_all.samples.rds")

#loaded_pos_table <- readRDS("path/to/your/directory/pos_table.rds")

##POSITIVE
A_B <- attr(pos_table, "intersections")$`A:B`
A_C <- attr(pos_table, "intersections")$`A:C`
A_D <- attr(pos_table, "intersections")$`A:D`
B_C <- attr(pos_table, "intersections")$`B:C`
B_D <- attr(pos_table, "intersections")$`B:D`
B_C_D <- attr(pos_table, "intersections")$`B:C:D`
A_B_C_D <- attr(pos_table, "intersections")$`A:B:C:D`



# Obtener la longitud máxima de los vectores
max_length <- max(length(A_B), length(A_C), length(A_D), length(B_C), length(B_D), length(B_C_D), length(A_B_C_D))

# Completar los vectores con NA para que tengan la misma longitud
A_B <- c(A_B, rep(0, max_length - length(A_B)))
A_C <- c(A_C, rep(0, max_length - length(A_C)))
A_D <- c(A_D, rep(0, max_length - length(A_D)))
B_C <- c(B_C, rep(0, max_length - length(B_C)))
B_D <- c(B_D, rep(0, max_length - length(B_D)))
B_C_D <- c(B_C_D, rep(0, max_length - length(B_C_D)))
A_B_C_D <- c(A_B_C_D, rep(0, max_length - length(A_B_C_D)))

# ZAP+
ZAP_POS <- data.frame(A_B, A_C, A_D, B_C, B_D, B_C_D, A_B_C_D)
names(ZAP_POS) <- c('SAMPLES1|2', 'SAMPLES1|3', 'SAMPLES1|4', 'SAMPLES2|3', 'SAMPLES2|4', 'SAMPLES2|3|4', 'ALL')

write.csv(ZAP_POS, file = "/home/adeluca/Documents/ZAP+_intersections.csv", row.names = F)

##NEGATIVE
A_B <- attr(neg_table, "intersections")$`A:B`
A_C <- attr(neg_table, "intersections")$`A:C`
A_D <- attr(neg_table, "intersections")$`A:D`
B_D <- attr(neg_table, "intersections")$`B:D`
C_D <- attr(neg_table, "intersections")$`C:D`

# Obtener la longitud máxima de los vectores
max_length <- max(length(A_B), length(A_C), length(A_D), length(B_D), length(C_D))

# Completar los vectores con NA para que tengan la misma longitud
A_B <- c(A_B, rep(0, max_length - length(A_B)))
A_C <- c(A_C, rep(0, max_length - length(A_C)))
A_D <- c(A_D, rep(0, max_length - length(A_D)))
B_D <- c(B_D, rep(0, max_length - length(B_D)))
C_D <- c(B_C_D, rep(0, max_length - length(B_C_D)))

# ZAP-
ZAP_NEG <- data.frame(A_B, A_C, A_D, B_D, C_D)
names(ZAP_NEG) <- c('SAMPLES1|2', 'SAMPLES1|3', 'SAMPLES1|4', 'SAMPLES2|4', 'SAMPLES3|4')
write.csv(ZAP_NEG, file = "/home/adeluca/Documents/ZAP-_intersections.csv", row.names = F)
'''

# ZAP+ = 'KLRC4', 'AC104078.1', 'NPDC1', 'MEX3D', 'AC097511.1', 'KRT73', 'LINC01980', 'LINC02611', 'HOPX', 'JAKMIP1', 'B3GAT1', 'ZAP70'

FeaturePlot(S1, features = c('KLRC4', 'AC104078.1', 'NPDC1', 'MEX3D', 'AC097511.1', 'KRT73', 'LINC01980', 'LINC02611', 'HOPX', 'JAKMIP1', 'B3GAT1', 'ZAP70'),order = T) # location in the clusters
FeaturePlot(S2, features = c('KLRC4', 'AC104078.1', 'NPDC1', 'MEX3D', 'AC097511.1', 'KRT73', 'LINC01980', 'LINC02611', 'HOPX', 'JAKMIP1', 'B3GAT1', 'ZAP70'),order = T) # location in the clusters
FeaturePlot(S3, features = c('KLRC4', 'AC104078.1', 'NPDC1', 'MEX3D', 'AC097511.1', 'KRT73', 'LINC01980', 'LINC02611', 'HOPX', 'JAKMIP1', 'B3GAT1', 'ZAP70'),order = T) # location in the clusters
FeaturePlot(S4, features = c('KLRC4', 'AC104078.1', 'NPDC1', 'MEX3D', 'AC097511.1', 'KRT73', 'LINC01980', 'LINC02611', 'HOPX', 'JAKMIP1', 'B3GAT1', 'ZAP70'),order = T) # location in the clusters

S1 <- readRDS("/home/adeluca/Documents/S1.rds")
S2 <- readRDS("/home/adeluca/Documents/S2.rds")
S3 <- readRDS("/home/adeluca/Documents/S3.rds")
S4 <- readRDS("/home/adeluca/Documents/S4.rds")

## Coexpression of ZAP70+ and ZAP70- with those DE-genes

#ZAP+ = 'KLRC4', 'NPDC1', 'MEX3D', 'KRT73', 'HOPX', 'JAKMIP1', 'B3GAT1', 'ZAP70', 'AC104078.1', 'LINC01980', 'LINC02611', 'AC097511.1'
FeaturePlot(S4, features = c("ZAP70", "KLRC4"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "NPDC1"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "MEX3D"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "KRT73"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "HOPX"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "JAKMIP1"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "B3GAT1"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "ZAP70"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "AC104078.1"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "LINC01980"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "LINC02611"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "AC097511.1"), blend = TRUE, order = T)

#ZAP- = 'PLBD1', 'VCAN', 'MNDA', 'CHIT1', 'SRGAP3', 'BTBD8', 'ZBTB20-AS1'
FeaturePlot(S4, features = c("ZAP70", "PLBD1"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "VCAN"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "MNDA"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "CHIT1"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "SRGAP3"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "BTBD8"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "ZBTB20-AS1"), blend = TRUE, order = T)

# INFORMATION
# IL2RB, CD8A, CD247, LAG3 and KLRK1 prognosis markers of CLL in paper
## file:///home/adeluca/Downloads/1471-2105-11-S9-S5.pdf

## Coexpression of ZAP70+ and ZAP70- with the common markers

#B-cells
FeaturePlot(S4, features = c("ZAP70", "MS4A1"), blend = TRUE, order = T)
#NK cells
FeaturePlot(S4, features = c("ZAP70", "GNLY"), blend = TRUE, order = T)

#T-celL
FeaturePlot(S4, features = c("ZAP70", "CD3E"), blend = TRUE, order = T)#
FeaturePlot(S4, features = c("ZAP70", "CD3G"), blend = TRUE, order = T)#
FeaturePlot(S4, features = c("ZAP70", "CD8A"), blend = TRUE, order = T)

#B-cells
FeaturePlot(S4, features = c("ZAP70", "CD79A"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "CD79B"), blend = TRUE, order = T)

#Monocytes
FeaturePlot(S4, features = c("ZAP70", "FCGR3A"), blend = TRUE, order = T)#

#Langerhans cells
FeaturePlot(S4, features = c("ZAP70", "LYZ"), blend = TRUE, order = T)

# Notably, ZAP-70 level in ALL is associated with CD38 expression, but no correlation was observed to specific cytogenetic abnormalities 
FeaturePlot(S4, features = c("ZAP70", "CD38"), blend = TRUE, order = T)


FeaturePlot(S4, features = c("ZAP70", "STAT4"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "THEMIS"), blend = TRUE, order = T)
FeaturePlot(S4, features = c("ZAP70", "HBA1"), blend = TRUE, order = T)

#ATM, BIRC3, EGR2, NFKBIE, NOTCH1, SF3B1, and TP53 putative drivers
##https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2023.1143811/full


## For instance, patients in subset #1 (utilize IGHV1/5/7 clan I genes, U-CLL) 
## and #2 (IGHV3-21/IGLV3-21, mixed SHM status) respond poorly to chemoimmunotherapy 
## and have a dismal outcome, whereas subset #4 patients (IGHV4-34/IGKV2-30, M-CLL) 
## show indolent disease courses and are rarely in need of treatment

## Coexpression of ZAP70+ and ZAP70- with Immunoglobulins

all.markers <- S1@assays$RNA@features@.Data
markers_names = row.names(all.markers)
AGUS = grep('CD38', markers_names)
ig_genes <- grep("^IG", markers_names, value = TRUE)

# IgVH - a reliable prognostic marker in CLL, equivalent to that of IgVH gene mutational status
# IgM - The expression levels of ZAP-70 in CLL cells are relatively stable over time. The aberrant ZAP-70 expression has further been found to associate with sIgM expression in CLL


## TAKE ALL MARKERS FROM ALL SAMPLES AND CLUSTERS WITH ZAP70 TO FIND A SIGNATURE OF MARKERS FOR ZAP70+ CELLS

mark1_ <- read.csv('/home/adeluca/Documents/S_180574/markers/markers_S180574.csv', header = F)
cluster4_S1 <- mark1_[mark1_$V7 == 4, ]
cluster4_S1 <- cbind(cluster4_S1, 'SAMPLE' = rep('1', 3736))
cluster2_S1 <- mark1_[mark1_$V7 == 2, ]
cluster2_S1 <- cbind(cluster2_S1, 'SAMPLE' = rep('1', 4771))

mark2_ <- read.csv('/home/adeluca/Documents/S_210183/markers/markers_S210183.csv', header = F)
cluster0_S2 <- mark2_[mark2_$V7 == 0, ]
cluster0_S2 <- cbind(cluster0_S2, 'SAMPLE' = rep('2', 2571))
cluster1_S2 <- mark2_[mark2_$V7 == 1, ]
cluster1_S2 <- cbind(cluster1_S2, 'SAMPLE' = rep('2', 2950))
cluster2_S2 <- mark2_[mark2_$V7 == 2, ]
cluster2_S2 <- cbind(cluster2_S2, 'SAMPLE' = rep('2', 5014))
cluster3_S2 <- mark2_[mark2_$V7 == 3, ]
cluster3_S2 <- cbind(cluster3_S2, 'SAMPLE' = rep('2', 377))
cluster4_S2 <- mark2_[mark2_$V7 == 4, ]
cluster4_S2 <- cbind(cluster4_S2, 'SAMPLE' = rep('2', 1748))
cluster5_S2 <- mark2_[mark2_$V7 == 5, ]
cluster5_S2 <- cbind(cluster5_S2, 'SAMPLE' = rep('2', 6486))

mark3_ <- read.csv('/home/adeluca/Documents/S_210414/markers/markers_S210414.csv', header = F)
cluster0_S3 <- mark3_[mark3_$V7 == 0, ]
cluster0_S3 <- cbind(cluster0_S3, 'SAMPLE' = rep('3', 543))
cluster1_S3 <- mark3_[mark3_$V7 == 1, ]
cluster1_S3 <- cbind(cluster1_S3, 'SAMPLE' = rep('3', 184))
cluster2_S3 <- mark3_[mark3_$V7 == 2, ]
cluster2_S3 <- cbind(cluster2_S3, 'SAMPLE' = rep('3', 2138))
cluster3_S3 <- mark3_[mark3_$V7 == 3, ]
cluster3_S3 <- cbind(cluster3_S3, 'SAMPLE' = rep('3', 10172))
cluster4_S3 <- mark3_[mark3_$V7 == 4, ]
cluster4_S3 <- cbind(cluster4_S3, 'SAMPLE' = rep('3', 749))

mark4_ <- read.csv('/home/adeluca/Documents/S_211855/markers/markers_S211855.csv', header = F)
cluster0_S4 <- mark4_[mark4_$V7 == 0, ]
cluster0_S4 <- cbind(cluster0_S4, 'SAMPLE' = rep('4', 2389))
cluster1_S4 <- mark4_[mark4_$V7 == 1, ]
cluster1_S4 <- cbind(cluster1_S4, 'SAMPLE' = rep('4', 3139))
cluster4_S4 <- mark4_[mark4_$V7 == 4, ]
cluster4_S4 <- cbind(cluster4_S4, 'SAMPLE' = rep('4', 3409))
cluster3_S4 <- mark4_[mark4_$V7 == 3, ]
cluster3_S4 <- cbind(cluster3_S4, 'SAMPLE' = rep('4', 7906))

all.markers <- rbind(cluster4_S1, cluster2_S1, cluster0_S2, cluster1_S2, cluster2_S2, cluster3_S2, cluster4_S2, cluster5_S2, cluster0_S3, cluster1_S3, cluster2_S3, cluster3_S3, cluster4_S3, cluster0_S4, cluster1_S4 , cluster4_S4, cluster3_S4)
names(all.markers) <- c("marker", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene", "sample")
all.markers <- all.markers[!duplicated(all.markers$marker), ]
write.csv(all.markers, file = "/home/adeluca/Documents/all.markers.csv", row.names = F)

all.markers <- read.csv('all.markers.csv')

top100_markers_ZAP70 <- head(all.markers, 500)


## Pearson correlation - gene to gene

seurat.data = S1[["RNA"]]$data
seurat.data.cor.pearson = cor(t(as.matrix(seurat.data)), 
                              method = "pearson")
# write.csv(seurat.data.cor, file =
# '../results/2019.12.16/E16.5_hipp_Seurat_correlations.csv'
# )

partial.coex.pearson = seurat.data.cor.pearson[rownames(seurat.data.cor.pearson) %in% 
                                                 c(tf1, tf2, hk), colnames(seurat.data.cor.pearson) %in% 
                                                 c(tf1, tf2, hk)]
diag(partial.coex.pearson) = 0


partial.coex.pearson = reshape2::melt(partial.coex.pearson)
colnames(partial.coex.pearson) = c("g1","g2","corr")

partial.coex.pearson$g1 <- factor(partial.coex.pearson$g1, c(tf1,hk,tf2))
partial.coex.pearson$g2 <- factor(partial.coex.pearson$g2, c(tf1,hk,tf2))

P = ggplot(partial.coex.pearson) + 
  geom_tile(aes(x=g1,y=g2, fill = corr),colour = "black", show.legend = TRUE) +
  #  facet_grid( g1 ~ g2  ,scales = "free", space = "free") + 
  scale_fill_gradient2(mid = "white",limits=c(-1, 1),low = "#DC0000B2", high = "#3C5488B2")+
  #scale_fill_gradient2(low = "darkred", mid = "white",  high = "darkblue", midpoint = 0,na.value = "grey80", space = "Lab", guide = "colourbar", aesthetics = "fill", limits = lim_coex, oob=scales::squish)+ theme(legend.position="bottom")+
  theme(#legend.title = element_blank(),
    #strip.text.x = element_text(color = "red"),
    #axis.text.y = element_text(color = ),
    axis.text.x = element_text(angle=45,hjust=1,vjust=1.0),
    legend.position="bottom"
  ) #+geom_text(aes(label=ifelse(t_hk == "hk", "H","")), color="grey", size=3)




correlation_matrix = cor(expr)