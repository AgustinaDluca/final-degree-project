# Load the object into R
data = readRDS('Documents/Documents/after_meeting/data.NOCD3/final/data_final.rds')

## Pipeline for RNA data
data@active.assay = 'RNA'
dim(data)
# 36601 20707

data$object <- 'SeuratObject'

# As we are not interested in keeping Tcells in my analysis we need to remove all those cells that express Tcells markers
# recognise CD3 cells 

CD3E <- names(which(data$RNA$counts.1['CD3E', ] >= 1)) # Get the cell names from the Seurat object
cell_names <- colnames(data) # Find the indices of the cells to remove
indices_to_remove <- which(cell_names %in% CD3E) # Remove the cells from the Seurat object
data <- data[, -indices_to_remove]

CD3E <- names(which(data$RNA$counts.2['CD3E', ] >= 1))
cell_names <- colnames(data)
indices_to_remove <- which(cell_names %in% CD3E)
data <- data[, -indices_to_remove]

CD3E <- names(which(data$RNA$counts.3['CD3E', ] >= 1))
cell_names <- colnames(data)
indices_to_remove <- which(cell_names %in% CD3E)
data <- data[, -indices_to_remove]

CD3E <- names(which(data$RNA$counts.4['CD3E', ] >= 1))
cell_names <- colnames(data)
indices_to_remove <- which(cell_names %in% CD3E)
data <- data[, -indices_to_remove]

#########################################################

CD3D <- names(which(data$RNA$counts.1['CD3D', ] >= 1))
cell_names <- colnames(data)
indices_to_remove <- which(cell_names %in% CD3D)
data <- data[, -indices_to_remove]

CD3D <- names(which(data$RNA$counts.2['CD3D', ] >= 1))
cell_names <- colnames(data)
indices_to_remove <- which(cell_names %in% CD3D)
data <- data[, -indices_to_remove]

CD3D <- names(which(data$RNA$counts.3['CD3D', ] >= 1))
cell_names <- colnames(data)
indices_to_remove <- which(cell_names %in% CD3D)
data <- data[, -indices_to_remove]

CD3D <- names(which(data$RNA$counts.4['CD3D', ] >= 1))
cell_names <- colnames(data)
indices_to_remove <- which(cell_names %in% CD3D)
data <- data[, -indices_to_remove]

#########################################################

CD3G <- names(which(data$RNA$counts.1['CD3G', ] >= 1))
cell_names <- colnames(data)
indices_to_remove <- which(cell_names %in% CD3G)
data <- data[, -indices_to_remove]

CD3G <- names(which(data$RNA$counts.2['CD3G', ] >= 1))
cell_names <- colnames(data)
indices_to_remove <- which(cell_names %in% CD3G)
data <- data[, -indices_to_remove]

CD3G <- names(which(data$RNA$counts.3['CD3G', ] >= 1))
cell_names <- colnames(data)
indices_to_remove <- which(cell_names %in% CD3G)
data <- data[, -indices_to_remove]

CD3G <- names(which(data$RNA$counts.4['CD3G', ] >= 1))
cell_names <- colnames(data)
indices_to_remove <- which(cell_names %in% CD3G)
data <- data[, -indices_to_remove]

## To observe how the removal went

dim(data)
# 36601 20356

VlnPlot(data, features = c('CD3E','CD3D', 'CD3G'), group.by = 'orig.ident')
VlnPlot(data, features = c('CD3E','CD3D', 'CD3G'), group.by = 'object')

Idents(data) = data$orig.ident
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

saveRDS(data, 'after_meeting/data.NOCD3/data_noCD3.rds')

## 2. NORMALIZATION

data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

## 3. IDENTIFY HIGHL VARIABLE FEATURES (feature selection)

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000) # 15-20% of the number of genes  TRY WITH 800

top10 <- head(VariableFeatures(data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## 4. SCALING - ensures that the magnitudes of gene expression are comparable across all genes, preventing genes with higher magnitudes from dominating the analysis.

data <- ScaleData(data, features = rownames(data)) # This function scales and centers the expression values for each gene across all cells. 

## 5. LINEAR DIMENSIONALITY REDUCTION (PCA)

data <- RunPCA(data, features = VariableFeatures(object = data))

print(data[['pca']], dims = 1:5, nfeatures = 5) # Examine and visualize PRINCIPAL COMPONENT results a few different ways

VizDimLoadings(data, dims = 1:2, reduction = 'pca') # shows the results of the PC's and the genes that are positive and negative. This function visualizes the loadings of features on the first two principal components. It helps you understand which genes contribute the most to the observed variation along these components.

DimPlot(data, reduction = 'pca') + NoLegend() # The position of cells along PC1 and PC2 represents their overall gene expression patterns with respect to the main sources of variation in the dataset. 
### SAVE ###
DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)
### SAVE ###
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)

## 6. DETERMINE THE DIMENSIONALITY OF THE dataSET

ElbowPlot(data)

## 7. CLUSTER THE CELLS

data <- FindNeighbors(data, dims = 1:20) #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.
data <- FindClusters(data, resolution = 0.5, cluster.name = 'clusters.res.0.5') #Partition the KNN graph into clusters.
table(data@active.ident)

#.   0    1    2    3    4    5    6    7    8     
#. 5509 5340 4801 1496 1438  845  644  204   79

## 8. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

data <- RunUMAP(data, dims = 1:20, reduction.name = 'umap.rna.res.0.5')

DimPlot(data, group.by = 'sample', reduction = 'umap.rna.res.0.5', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))
VlnPlot(data[, data$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'sample', log = TRUE, sort = T, cols = c('mediumseagreen','darkseagreen','seagreen','mediumaquamarine'))

DimPlot(data, group.by = 'matchScore_RT_annotation', reduction = 'umap.harmony.no.lymph6', label = F, cols = c('gray','gray','gray','seagreen','royalblue', 'gray', 'red', 'yellow', 'white','gray'))
VlnPlot(data[, data$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'matchScore_RT_annotation', log = TRUE, sort = T, cols = c('gray','red','yellow','seagreen','gray', 'royalblue', 'gray', 'white','gray'))

DimPlot(data, group.by = 'clusters.res.0.5', reduction = 'umap.rna.res.0.5', label = T, cols = c('mediumaquamarine','seagreen','mediumseagreen','darkseagreen','lightblue', 'cornflowerblue', 'royalblue', 'pink', 'palevioletred'))
VlnPlot(data[, data$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'clusters.res.0.5', log = TRUE, sort = T, cols = c('mediumaquamarine','seagreen','darkseagreen','royalblue', 'lightblue', 'pink', 'palevioletred', 'cornflowerblue','gray'))

p1 + p2 

proportions <- as.data.frame((table(data$clusters.res.0.5, data$ZAP_expression) / rowSums(table(data$clusters.res.0.5, data$ZAP_expression)))*100)
# Assuming 'proportions' is your S1 frame
ZAP70_expression_plot <- ggplot(proportions, aes(x = factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Proportion of ZAP+ and ZAP- in each cluster",
       x = "Clusters",
       y = "Percentage",
       fill = "Expression") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  scale_fill_manual(values = c("ZAP-" = "skyblue", "ZAP+" = "salmon"))

ZAP70_expression_plot

VlnPlot(data, features = 'nFeature_RNA', group.by = 'clusters.res.0.5', log = TRUE, sort = T, cols = c('palevioletred', 'pink', 'lightblue', 'cornflowerblue', 'royalblue', 'darkseagreen', 'seagreen', 'mediumaquamarine', 'mediumseagreen'))

FeaturePlot(data, features = 'ZAP70', order  = TRUE, cols = c("grey", "red"), pt.size = 0.5, reduction = 'umap.rna.res.0.5')

# INTEGRATIONS (https://satijalab.org/seurat/articles/seurat5_integration)
data = SetIdent(data, value = data$orig.ident)

# 1. Anchor-based CCA integration (method=CCAIntegration)
data <- IntegrateLayers(object = data, method = CCAIntegration,
                        orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

data <- FindNeighbors(data, reduction = "integrated.cca", dims = 1:30)
data <- FindClusters(data, resolution = 0.5, cluster.name = "cca_clusters")
data <- RunUMAP(data, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
DimPlot(data,reduction = "umap.cca",group.by = "cca_clusters", )

DimPlot(data, reduction = "umap.cca", group.by = "clusters.res.0.5", , cols = c('mediumaquamarine','seagreen','mediumseagreen','darkseagreen','lightblue', 'cornflowerblue', 'royalblue', 'pink', 'palevioletred'))

# Samples
DimPlot(data, reduction = "umap.cca", group.by = "orig.ident",  cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))
# RT annotation
DimPlot(data, reduction = "umap.cca", group.by = "matchScore_RT_annotation",  cols = c('gray', 'gray', 'gray', 'seagreen', 'royalblue', 'gray', 'red', 'yellow', 'white', 'gray'))
# ZAP70 expression
FeaturePlot(data, features = 'ZAP70', reduction = 'umap.cca')


# 2. Anchor-based RPCA integration (method=RPCAIntegration)
data <- IntegrateLayers(object = data, method = RPCAIntegration,
                        orig.reduction = "pca", new.reduction = "integrated.rpca",
                        verbose = FALSE)

data <- FindNeighbors(data, reduction = "integrated.rpca", dims = 1:30)
data <- FindClusters(data, resolution = 0.5, cluster.name = "rpca_clusters")
data <- RunUMAP(data, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
DimPlot(data,reduction = "umap.rpca",group.by = "rpca_clusters",combine = FALSE)

DimPlot(data, reduction = "umap.rpca", group.by = "clusters.res.0.5",  cols = c('mediumaquamarine','seagreen','mediumseagreen','darkseagreen','lightblue', 'cornflowerblue', 'royalblue', 'pink', 'palevioletred'))

# Samples
DimPlot(data, reduction = "umap.rpca", group.by = "orig.ident",  cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))
# RT annotation
DimPlot(data, reduction = "umap.rpca", group.by = "matchScore_RT_annotation",  cols = c('gray', 'gray', 'gray', 'seagreen', 'royalblue', 'gray', 'red', 'yellow', 'black', 'gray'))
# ZAP70 expression
FeaturePlot(data, features = 'ZAP70', reduction = 'umap.rpca')

# 3. Harmony (method=HarmonyIntegration)
data <- IntegrateLayers(object = data, method = HarmonyIntegration,
                        orig.reduction = "pca", new.reduction = "harmonyoriginal",
                        verbose = FALSE)

data <- FindNeighbors(data, reduction = "harmonyoriginal", dims = 1:30)
data <- FindClusters(data, resolution = 0.5, cluster.name = "harmony_original_clusters")
data <- RunUMAP(data, reduction = "harmonyoriginal", dims = 1:30, reduction.name = "umap.harmony.original")
DimPlot(data,reduction = "umap.harmony.original",group.by = "harmony_original_clusters",combine = FALSE)

DimPlot(data, reduction = "umap.harmony.original", group.by = "clusters.res.0.5",  cols = c('mediumaquamarine','seagreen','mediumseagreen','darkseagreen','lightblue', 'cornflowerblue', 'royalblue', 'pink', 'palevioletred'))

# Samples
DimPlot(data, reduction = "umap.harmony.original", group.by = "orig.ident",  cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))
# RT annotation
DimPlot(data, reduction = "umap.harmony.original", group.by = "matchScore_RT_annotation",  cols = c('gray', 'gray', 'gray', 'seagreen', 'royalblue', 'gray', 'red', 'yellow', 'black', 'gray'))
# ZAP70 expression
FeaturePlot(data, features = 'ZAP70', reduction = 'umap.harmony.original')

# 4.FastMNN (method= FastMNNIntegration)
data <- IntegrateLayers(object = data, method = FastMNNIntegration,
                        new.reduction = "integrated.mnn",
                        verbose = FALSE)

data <- FindNeighbors(data, reduction = "integrated.mnn", dims = 1:30)
data <- FindClusters(data, resolution = 0.5, cluster.name = "mnn_clusters")
data <- RunUMAP(data, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")
DimPlot(data,reduction = "umap.mnn",group.by = "mnn_clusters",combine = FALSE)

DimPlot(data, reduction = "umap.mnn", group.by = "clusters.res.0.5",  cols = c('mediumaquamarine','seagreen','mediumseagreen','darkseagreen','lightblue', 'cornflowerblue', 'royalblue', 'pink', 'palevioletred'))

# Samples
DimPlot(data, reduction = "umap.mnn", group.by = "orig.ident",  cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))
# RT annotation
DimPlot(data, reduction = "umap.mnn", group.by = "matchScore_RT_annotation",  cols = c('gray', 'gray', 'gray', 'seagreen', 'royalblue', 'gray', 'red', 'yellow', 'black', 'gray'))
# ZAP70 expression
FeaturePlot(data, features = 'ZAP70', reduction = 'umap.mnn')

## RUN HARMONY NORMAL
data <- RunHarmony(data, "sample")
data <- RunUMAP(data, dims = 1:20, reduction = 'harmony', reduction.name = 'umap.harmony')

data <- FindNeighbors(data, reduction = "harmony", dims = 1:30)
data <- FindClusters(data, resolution = 0.5, cluster.name = "harmony_clusters")

DimPlot(data, reduction = "umap.harmony", group.by = "clusters.res.0.5",  cols = c('mediumaquamarine','seagreen','mediumseagreen','darkseagreen','lightblue', 'cornflowerblue', 'royalblue', 'pink', 'palevioletred'))
DimPlot(data, reduction = "umap.harmony", group.by = "matchScore_RT_annotation",  cols = c('gray','gray','gray','seagreen','royalblue', 'gray', 'red', 'yellow', 'black', 'gray'))
DimPlot(data,reduction = "umap.harmony",group.by = "sample", cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))
DimPlot(data,reduction = "umap.harmony",group.by = "harmony_clusters")


FeaturePlot(data, features = 'ZAP70', reduction = 'umap.harmony')

## barplot to observe the proportion of RT_annotations per cluster
proportions = as.data.frame((table(data$harmony_clusters, data$matchScore_RT_annotation) / rowSums(table(data$harmony_clusters, data$matchScore_RT_annotation)))*100) 

ggplot(proportions, aes(x = as.factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +  # Change position to "fill"
  labs(x = "Clusters", y = "Proportion", fill = "matchScore_RT_annotation") +
  theme_bw() +
  theme(legend.position = "bottom") +   
  scale_fill_manual(values = c("lightblue","cornflowerblue", "blue4","seagreen","royalblue","palegoldenrod", "red","yellow",'black', "gray"))

#########################################################

## A. TRYING OTHER RESOLUTIONS (0.8)

data <- FindClusters(data, resolution = 0.8) #Partition the KNN graph into clusters.
table(data@active.ident)
#   0    1    2    3    4    5    6    7    8    9   
# 4868 4261 3361 3044 2064 1850  437  267  104  100 

data <- RunUMAP(data, dims = 1:20, reduction.name = 'umap.rna.0.8')

p1 = DimPlot(data, group.by = 'RNA_snn_res.0.8', reduction = 'umap.harmony', label = T)
p2 = DimPlot(data, group.by = 'matchScore_RT_annotation', reduction = 'umap.harmony', label = F, cols = c('gray','gray','gray','seagreen','royalblue', 'gray', 'red', 'yellow', 'black','gray'))

p1 | p2

## B. TRYING OTHER RESOLUTIONS (1)

data <- FindClusters(data, resolution = 1) #Partition the KNN graph into clusters.
table(data@active.ident)

#   0    1    2    3    4    5    6    7    8    9   10   11  
# 4247 3857 3431 2268 2055 1693 1528  441  356  267  109  104 

data <- RunUMAP(data, dims = 1:20, reduction.name = 'umap.rna.1')

p1 = DimPlot(data, group.by = 'RNA_snn_res.1', reduction = 'umap.harmony', label = T)
p2 = DimPlot(data, group.by = 'matchScore_RT_annotation', reduction = 'umap.harmony', label = F, cols = c('gray','gray','gray','seagreen','royalblue', 'gray', 'red', 'yellow', 'black','gray'))
p3 = DimPlot(data, group.by = 'orig.ident', reduction = 'umap.harmony', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))

p1 | p2 | p3

## B. TRYING OTHER RESOLUTIONS (1.2)

data <- FindClusters(data, resolution = 1.2) #Partition the KNN graph into clusters.
table(data@active.ident)

#   0    1    2    3    4    5    6    7    8    9   10   11   12    
# 4148 3492 3347 2368 2041 2006 1576  443  315  277  128  111  104

data <- RunUMAP(data, dims = 1:20, reduction.name = 'umap.rna.1.2')

p1 = DimPlot(data, group.by = 'RNA_snn_res.1.2', reduction = 'umap.harmony', label = T)
p2 = DimPlot(data, group.by = 'matchScore_RT_annotation', reduction = 'umap.harmony', label = F, cols = c('gray','gray','gray','seagreen','royalblue', 'gray', 'red', 'yellow', 'black','gray'))
p3 = DimPlot(data, group.by = 'orig.ident', reduction = 'umap.harmony', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))

p1 | p2 | p3

#####################
# Find the best resolution value and integration method

p1 = DimPlot(data, group.by = 'RNA_snn_res.1', reduction = 'umap.harmony', label = T)
# matchscore
p2 = DimPlot(data, group.by = 'matchScore_RT_annotation', reduction = 'umap.harmony', label = F, cols = c('gray','gray','gray','seagreen','royalblue', 'gray', 'red', 'yellow', 'black','gray'))
# samples
p3 = DimPlot(data, group.by = 'orig.ident', reduction = 'umap.harmony', label = T)

proportions = as.data.frame((table(data$RNA_snn_res.1, data$matchScore_RT_annotation) / rowSums(table(data$RNA_snn_res.1, data$matchScore_RT_annotation)))*100) 
p4 = ggplot(proportions, aes(x = as.factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +  # Change position to "fill"
  labs(x = "Clusters", y = "Proportion", fill = "matchScore_RT_annotation") +
  theme_bw() +
  theme(legend.position = "bottom") +   
  scale_fill_manual(values = c("gray","gray", "gray","seagreen","royalblue", "gray", "red","yellow",'black', "gray"))

p1 | p4 | p2


# unify clusters nearby
table(data$RNA_snn_res.1)

library(plyr)
data$my_clusters <- revalue(data$RNA_snn_res.1, c("0"="0", "1"="1", "10"= "1", "11"= "7", "2" = "0","3" = "2","4" = "3","5" = "2","6" = "4","7" = "5", "8" = "4","9" = "6"))

p1 = DimPlot(data, group.by = 'my_clusters', reduction = 'umap.harmony', label = T)
# matchscore
p2 = DimPlot(data, group.by = 'matchScore_RT_annotation', reduction = 'umap.harmony', label = F, cols = c('gray','gray','gray','seagreen','royalblue', 'gray', 'red', 'yellow', 'black','gray'))

proportions = as.data.frame((table(data$my_clusters, data$matchScore_RT_annotation) / rowSums(table(data$my_clusters, data$matchScore_RT_annotation)))*100) 
p4 = ggplot(proportions, aes(x = as.factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +  # Change position to "fill"
  labs(x = "Clusters", y = "Proportion", fill = "matchScore_RT_annotation") +
  theme_bw() +
  theme(legend.position = "bottom") +   
  scale_fill_manual(values = c("gray","gray", "gray","seagreen","royalblue", "gray", "red","yellow",'black', "gray"))

p1 | p4 | p2

# LOOK FOR TCELL MARKERS
p5 = DotPlot(data, features = c('CXCR4', 'CD24','CD27', 'MIR155HG','CCND2', 'PCNA', 'MKI67','MZB1','IGHM', 'XBP1'), group.by = 'my_clusters') + RotatedAxis()

p1 | p5

FeaturePlot(data, features = c('CD3E', 'CD3G', 'CD8A', 'IL7R'), reduction = 'umap.harmony')

#####################
library(reshape2)
## By sample, prognosis, ZAP70 expression and clusters
proportions <- as.data.frame((table(data$orig.ident, data$ZAP_expression, data$prognosis, data$new_annotation) / rowSums(table(data$orig.ident, data$ZAP_expression, data$prognosis, data$new_annotation)))*100)
long_proportions <- melt(proportions, varnames = c("Sample","ZAP_expression", "prognosis", "raw_annotation_clusters"), value.name = "proportion")

ggplot(long_proportions, aes(x = Var4, y = proportion, fill = interaction(Var3, Var2))) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Var1, scales = "free_x") +
  labs(x = "Cell type raw annotation", y = "Proportion (%)", fill = "Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  scale_fill_manual(values = cols_cond)


library(ggalluvial)

# Convert necessary columns to factors for better plotting
long_proportions$Var1 <- as.factor(long_proportions$Var1) # sample
long_proportions$Var2 <- as.factor(long_proportions$Var2) # ZAP expression
long_proportions$Var3 <- as.factor(long_proportions$Var3) # prognosis
long_proportions$Var4 <- as.factor(long_proportions$Var4) # raw_annotation_clusters

# Create the alluvial plot
ggplot(long_proportions,
       aes(axis1 = Var4, axis2 = Var3, axis3 = Var2, y = proportion)) +
  geom_alluvium(aes(fill = Var4), width = 1/12) + # Change to Var4 for clusters
  geom_stratum(width = 1/12, fill = "lightgray", color = "gray") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Cell type", "Prognosis", "ZAP70 Expression"), expand = c(0.15, 0.05)) +
  scale_fill_manual(values = cols_anno) + # Ensure cols_cond has cluster colors
  labs(title = "Alluvial Plot of ZAP70 Expression, Prognosis, and Raw Annotation Clusters",
       x = "Categories", y = "Proportion (%)", fill = "Clusters") + # Change legend title to "Clusters"
  theme_minimal() +
  theme(legend.position = "top")

###############################
#I did this to JoinLayers() in order to calculate the markers and see how the removal of CD3 cells affected and also retest the model of the RTannotation
data1 = data
data1 = SetIdent(data1, value = data1$orig.ident)
data1 = JoinLayers(data1)

Idents(data1) = data1$my_clusters

markers <- FindAllMarkers(data1, only.pos = T)
write.csv(markers, file = 'after_meeting/data.NOCD3/markers_myclusters.csv')
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

DoHeatmap(data1, features = top10$gene) + NoLegend()

top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=100)
write.csv(top, file = 'after_meeting/data.NOCD3/top100markers_myclusters.csv')

top_table = as.data.frame(top)

top_table$X0
# I consider important: IGKC, GNG7, JCHAIN, MAST4, CCND3, NOTCH2 (cluster 'MIR155HGhi' MIR155HG)

top_table$X1
# I consider important: IGLC2, IGLC3, IGHM, CXCR4, JUND, JUNB, CD74, CCDC50, CCDC50 (RPS...genes)

top_table$X2
# I consider important: MIR3681HG, CCDC26 (AC...genes)

top_table$X3
# I consider important: CD69, CCDC50, XBP1, RELB

top_table$X4
# I consider important: NOTCH2, GNG7, MAST4, JCHAIN (AC...genes)

top_table$X5
# I consider important: CD69, CBX4, ZAP70, DNAJB11, XBP1

top_table$X6
# I consider important: TNFAIP3, STAT4, IL7R, GZMH, ZAP70, CD28, THEMIS, GNLY, NKG7, CD8A, FOS (Lymohid cluster?)


data@active.assay = 'RNA'
DotPlot(data, features = c('CXCR4', 'CD24','CD27', 'MIR155HG','CCND2', 'PCNA', 'MKI67','MZB1','IGHM', 'XBP1'), group.by = 'raw_annotation_clusters', dot.scale = 8) + RotatedAxis()


model = readRDS('after_meeting/RT_model/model_RT_annotations.rds')
data_rt = readRDS('RT/data_rt.rds')

Idents(data_rt) = data_rt$annotation_final
markers <- FindAllMarkers(data_rt, only.pos = T)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=100)

top = as.data.frame(top)
## Cell projection in my data
out <- identity_map(scale.data = data1@assays$RNA$scale.data,model = model,top)

### cell identities
ids <- out$ids 

data1$agustina = ids

p1 = DimPlot(data, reduction = "umap.harmony", group.by = "matchScore_RT_annotation", cols = c('gray','gray','gray','seagreen','palegreen3', 'gray', 'red', 'yellow', 'black', 'gray'))
p2 = DimPlot(data1, reduction = "umap.harmony", group.by = "agustina", cols = c('gray','gray','gray','seagreen','palegreen3', 'gray', 'red', 'yellow', 'black', 'gray'))

################
# Remove cluster 6 as it is a lymphoid cluster
cluster6 <- which(data@meta.data$my_clusters %in% '6') # identify T cells (positions)
length(cluster6) # 267

data <- subset(data, cells = -cluster6) # subset to remove T cells

## CLUSTER THE CELLS
data <- FindNeighbors(data, dims = 1:20) #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.
data <- FindClusters(data, resolution = 0.5, cluster.name = 'clusters.no.lymph6') #Partition the KNN graph into clusters.
table(data@active.ident)

#.   0    1    2    3    4    5    6    7        
#. 5506 5179 4844 1454 1436  995  596   79

## 8. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

data <- RunUMAP(data, dims = 1:20, reduction.name = 'umap.no.lymph6')
DimPlot(data, group.by = 'clusters.no.lymph6', reduction = 'umap.no.lymph6', label = T)

## RUN HARMONY NORMAL
data <- RunHarmony(data, "sample")
data <- RunUMAP(data, dims = 1:20, reduction = 'harmony', reduction.name = 'umap.harmony.no.lymph6')

data <- FindNeighbors(data, reduction = "harmony", dims = 1:30)
data <- FindClusters(data, resolution = 0.5, cluster.name = "harmony.no.lymph6")

DimPlot(data, reduction = "umap.harmony.no.lymph6", group.by = "clusters.no.lymph6",  cols = c('mediumaquamarine','seagreen','mediumseagreen','darkseagreen','lightblue', 'cornflowerblue', 'royalblue', 'pink', 'palevioletred'))
DimPlot(data, reduction = "umap.harmony.no.lymph6", group.by = "matchScore_RT_annotation",  cols = c('gray','gray','gray','seagreen','royalblue', 'gray', 'red', 'yellow', 'black', 'gray'))
DimPlot(data,reduction = "umap.harmony.no.lymph6",group.by = "sample", cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))
DimPlot(data,reduction = "umap.harmony.no.lymph6",group.by = "harmony.no.lymph6")

FeaturePlot(data, features = 'ZAP70', reduction = 'umap.harmony.no.lymph6')

## Increase the resolution up to 1

data <- FindClusters(data, resolution = 1) #Partition the KNN graph into clusters.
table(data@active.ident)

#   0    1    2    3    4    5    6    7    8    9   10   11  
# 4247 3857 3431 2268 2055 1693 1528  441  356  267  109  104 

data <- RunUMAP(data, dims = 1:20, reduction.name = 'umap.rna.1.no6')

p1 = DimPlot(data, group.by = 'RNA_snn_res.1', reduction = 'umap.harmony.no.lymph6', label = T)
p2 = DimPlot(data, group.by = 'matchScore_RT_annotation', reduction = 'umap.harmony.no.lymph6', label = F, cols = c('gray','gray','gray','seagreen','royalblue', 'gray', 'red', 'yellow', 'black','gray'))
p3 = DimPlot(data, group.by = 'orig.ident', reduction = 'umap.harmony.no.lymph6', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))

p1 | p2 | p3


#####################

# unify clusters nearby
table(data$RNA_snn_res.1)

library(plyr)
data$my_final_clusters <- revalue(data$RNA_snn_res.1, c("0"="0", "1"="1", "10"= "1", "11"= "2", "2" = "0","3" = "2","4" = "3","5" = "2","6" = "4","7" = "5", "8" = "4","9" = "6"))
data$my_final_clusters <- revalue(data$my_final_clusters, c("0"="0", "1"="1", "2"= "1", "3"= "2", "4" = "3","5" = "4","6" = "5"))

p1 = DimPlot(data, group.by = 'my_final_clusters', reduction = 'umap.harmony.no.lymph6', label = F, cols = c('#8bc28c', '#f88f58','#f18aad', '#ea6759', '#6667ab', '#f3c65f'))
# matchscore
p2 = DimPlot(data, group.by = 'matchScore_RT_annotation', reduction = 'umap.harmony.no.lymph6', label = F, cols = c("gray","gray", "gray","#297d4e","#6d8ce8", "gray", "#ff4c4c","#ffff66",'#46494c', "gray"))

proportions = as.data.frame((table(data$my_final_clusters, data$matchScore_RT_annotation) / rowSums(table(data$my_final_clusters, data$matchScore_RT_annotation)))*100) 
p4 = ggplot(proportions, aes(x = as.factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +  # Change position to "fill"
  labs(x = "Clusters", y = "Proportion", fill = "matchScore_RT_annotation") +
  theme_bw() +
  theme(legend.position = "bottom") +   
  scale_fill_manual(values = c("gray","gray", "gray","#297d4e","#6d8ce8", "gray", "#ff4c4c","#ffff66",'#46494c', "gray"))

p1 | p4 | p2

## ZAP70

FeaturePlot(data, features = 'ZAP70', reduction = 'umap.harmony.no.lymph6')
DimPlot(data, group.by = 'clusters.res.0.5', reduction = 'umap.harmony.no.lymph6', label = T, cols = c('mediumaquamarine','seagreen','mediumseagreen','darkseagreen','lightblue', 'cornflowerblue', 'royalblue', 'pink', 'palevioletred'))
DimPlot(data, group.by = 'orig.ident', reduction = 'umap.harmony.no.lymph6', label = T, cols = c('mediumaquamarine','seagreen','mediumseagreen','darkseagreen','lightblue', 'cornflowerblue', 'royalblue', 'pink', 'palevioletred'))

# LOOK FOR TCELL MARKERS
p5 = DotPlot(data, features = c('CXCR4', 'CD24','CD27', 'MIR155HG','CCND2', 'PCNA', 'MKI67','MZB1','IGHM', 'XBP1'), group.by = 'my_final_clusters') + RotatedAxis()

p1 | p5

# COLOR CELLS THAT ARE ZAP+ 
DimPlot(data, group.by = "ZAP_expression", cols = c("lightgray", "red"), pt.size = 0.1, order = T, reduction = 'umap.harmony.no.lymph6')
table(data$ZAP_expression)

# EXPRESSION
FeaturePlot(data, features = 'ZAP70', order  = TRUE, cols = c("grey", "red"), pt.size = 0.5, reduction = 'umap.harmony.no.lymph6')

VlnPlot(data, group.by = 'orig.ident', features = 'ZAP70', sort = T)
VlnPlot(data, features = 'ZAP70', group.by = 'my_final_clusters',slot = "counts", log = TRUE) # RNA counts

proportions <- as.data.frame((table(data$my_final_clusters, data$ZAP_expression) / rowSums(table(data$my_final_clusters, data$ZAP_expression)))*100)
ZAP70_expression_plot <- ggplot(proportions, aes(x = factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Proportion of ZAP+ and ZAP- in each cluster",
       x = "Clusters",
       y = "Percentage",
       fill = "Expression") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  scale_fill_manual(values = c("ZAP-" = "skyblue", "ZAP+" = "salmon"))

ZAP70_expression_plot

## RELATED MARKERS w/ ZAP70 from the RT dataset

FeaturePlot(data, features = c('CXCR4', 'CD27', 'IGHM'), order = T, max.cutoff = 'q90', label = F, reduction = 'umap.harmony.no.lymph6') # location in the clusters

p1 = FeaturePlot(data, features = 'CXCR4', order = T, max.cutoff = 'q90', label = T, reduction = 'umap.harmony.no.lymph6') # location in the clusters
p2 = VlnPlot(data, features = 'CXCR4', sort = T)
p1*p2

p1 = FeaturePlot(data, features = 'CD27', order = T, max.cutoff = 'q90', label = T, reduction = 'umap.harmony.no.lymph6') # location in the clusters
p2 = VlnPlot(data, features = 'CD27', sort = T)
p1*p2

p1 = FeaturePlot(data, features = 'IGHM', order = T, max.cutoff = 'q90', label = T, reduction = 'umap.harmony.no.lymph6') # location in the clusters
p2 = VlnPlot(data, features = 'IGHM', sort = T)
p1*p2

# Coexpression with ZAP70
FeaturePlot(data, features = c("ZAP70", "CXCR4"), blend = TRUE, order = T, reduction = 'umap.harmony.no.lymph6')
FeaturePlot(data, features = c("ZAP70", "CD27"), blend = TRUE, order = T, reduction = 'umap.harmony.no.lymph6')
FeaturePlot(data, features = c("ZAP70", "IGHM"), blend = TRUE, order = T, reduction = 'umap.harmony.no.lymph6')

VlnPlot(data, features = c('CXCR4', 'ZAP70', 'CD27', 'IGHM'), group.by = 'my_final_clusters', sort = T)

VlnPlot(data, features = c('CXCR4', 'CD27'), group.by = 'ZAP_expression')

VlnPlot(data, features = c('CXCR4', 'CD27', 'MIR155HG', 'TP53INP1', 'CCND2'), group.by = 'seurat_clusters', split.by = 'ZAP_expression')
# rosa ZAP70-
# azul ZAP70+

#####################

#I did this to JoinLayers() in order to calculate the markers and see how the removal of CD3 cells affected and also retest the model of the RTannotation
data1 = data
data1 = SetIdent(data1, value = data1$orig.ident)
data1 = JoinLayers(data1)

Idents(data1) = data1$my_final_clusters

markers <- FindAllMarkers(data, only.pos = T)
write.csv(markers, file = 'after_meeting/data.NOCD3/markers_myclusters.csv')
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

DoHeatmap(data1, features = top10$gene) + NoLegend()

top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=100)
write.csv(top, file = 'after_meeting/data.NOCD3/top100markers_myclusters.csv')

top_table = as.data.frame(top)

top_table$X0
# I consider important: IGKC, GNG7, JCHAIN, MAST4, CCND3, NOTCH2 

top_table$X1
# I consider important: IGLC2, IGLC3, IGHM, CXCR4, JUND, JUNB, CCDC50, HLA-B, XBP1

top_table$X2
# I consider important: IGLC2, IGHM, IGLC3, MIR3681HG, JUND, CD69 (RPS...genes)

top_table$X3
# I consider important: CD69, CCDC50, XBP1, RELB, IGLC3, FOSB, CD69, JUND, 

top_table$X4
# I consider important: MAST4, NOTCH2, GNG7, MAST4, JCHAIN (AC...genes)

top_table$X5
# I consider important: ZAP70, DNAJB9, CD69, DNAJB11, XBP1

top_table$X6
# I consider important: 

DotPlot(data1, features = c('CXCR4', 'CD24','CD27', 'MIR155HG','CCND2', 'PCNA', 'MKI67','MZB1','IGHM', 'XBP1', 'IGLC3'), group.by = 'my_final_clusters', col.min = 0) + RotatedAxis()
####################

# Cluster annotation (my_final_clusters)
# with the help of the rtdatasetmarkers, myfinalclustersmarkers, bibliography I annotated the clusters

data$my_final_clusters <- revalue(data$my_final_clusters, c("0"="0", "1"="1", "2"= "1", "3"= "2", "4" = "3","5" = "4","6" = "5"))

Idents(data) = data$my_final_clusters
markers <- FindAllMarkers(data, only.pos = T)
write.csv(markers, file = 'Documents/Documents/after_meeting/data.NOCD3/final/markers_myfinalclusters.csv')

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=100)

features <- c('CXCR4','CD24', 'CD27', 'TXNRD1' , 'MAST4', 'NOTCH2', 'MZB1', 'IGHM','XBP1', 'GTF2E2', 'WDFY1', 'HDAC9')
DotPlot(data, features = features, group.by = 'new_annotation', dot.scale = 8) +
  scale_color_gradient(low = "yellow", high = "purple") +
  RotatedAxis()

VlnPlot(data, features = 'ZAP70', sort = T, log = T, slot = 'scale.data', cols = cols_anno)

DotPlot(data, features = 'ZAP70', group.by = 'new_annotation', dot.scale = 8) +
  scale_color_gradient(low = "yellow", high = "purple") +
  RotatedAxis() + coord_flip()

proportions <- as.data.frame((table(data$new_annotation, data$ZAP_expression) / rowSums(table(data$new_annotation, data$ZAP_expression)))*100)
# Assuming 'proportions' is your S1 frame
ZAP70_expression_plot <- ggplot(proportions, aes(x = factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Proportion of ZAP+ and ZAP- in each cluster",
       x = "Clusters",
       y = "Percentage",
       fill = "Expression") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  scale_fill_manual(values = c("ZAP-" = "#eee3ab", "ZAP+" = "#80a1c1")) + coord_flip()

ZAP70_expression_plot


## HEATMAP
features <- c('CXCR4','CD24', 'CD27', 'TXNRD1' , 'MAST4', 'NOTCH2', 'MZB1', 'IGHM','XBP1' )

data_h <- FetchData(data, vars = features, slot = "scale.data")
cluster_info <- data@meta.data$new_annotation

# Aggregate data by clusters
data_aggregated <- aggregate(. ~ cluster_info, data_h, mean)
rownames(data_aggregated) <- data_aggregated$cluster_info
data_matrix <- as.matrix(data_aggregated[, -1])  # Remove the cluster column

# Get the order of rows based on cluster labels
cluster_order <- order(data_aggregated$cluster_info, decreasing = T)

# Reorder the rows of data_matrix
data_matrix_ordered <- data_matrix[cluster_order, ]

# Plot the heatmap with ordered rows
pheatmap(data_matrix_ordered,
         cluster_rows = FALSE,  # Since rows are already ordered
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         fontsize_col = 10,
         angle_col = 45, scale = 'column')

data <- SetIdent(data, value = data$my_final_clusters)
new_annotation <- c("CXCR4loIGHMlo", "CXCR4loTXNRD1hi", "CXCR4loIGHMmidMAST4hi","CXCR4midIGHMlo", "CXCR4midIGHMloNOTCH2mid", "CXCR4hiIGHMhi|CD27hiMZB1hiXBP1hi")

new_annotation <- new_annotation[as.factor(levels(data))]
data$new_annotation <- new_annotation[as.factor(data$my_final_clusters)]

names(new_annotation) <- levels(data)

data <- SetIdent(data, value = data$new_annotation)

DimPlot(data, reduction = 'umap.harmony.no.lymph6', label = F, cols = c('orange','palegreen3','seagreen','palegoldenrod','pink3', 'pink', 'palevioletred'))

markers <- FindAllMarkers(my_data, only.pos = T)
write.csv(markers, file = 'Documents/Documents/Supplementary material/RNA/markers/markers_new_annotation.csv')
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

DoHeatmap(my_data, features = top10$gene) + NoLegend()

top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=100)
write.csv(top, file = 'Documents/Documents/Supplementary material/RNA/markers/markers_new_annotation.csv')

### Enrichment GO pathways

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
## CLUSTERS (RNA)

CLL_transition <- ClosestFeature(data, regions = top$`CLL transition`)
# Perform pathway enrichment
enrich_result_CLL_transition <- enrichGO(gene = CLL_transition$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.5, qvalueCutoff = 0.5)
# Visualize the enriched pathways
dotplot(enrich_result_CLL_transition)
barplot(enrich_result_CLL_transition)

CXCR4hiCD27hiIGMLhi <- ClosestFeature(data, regions = top$`CXCR4hiCD27hiIGM|Lhi`)
# Perform pathway enrichment
enrich_result_CXCR4hiCD27hiIGMLhi <- enrichGO(gene = CXCR4hiCD27hiIGMLhi$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
# Visualize the enriched pathways
dotplot(enrich_result_CXCR4hiCD27hiIGMLhi)
barplot(enrich_result_CXCR4hiCD27hiIGMLhi)

CXCR4hiCD27midIGMLhi <- ClosestFeature(data, regions = top$`CXCR4hiCD27midIGM|Lhi`)
# Perform pathway enrichment
enrich_result_CXCR4hiCD27midIGMLhi <- enrichGO(gene = CXCR4hiCD27midIGMLhi$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.4, qvalueCutoff = 0.4)
# Visualize the enriched pathways
dotplot(enrich_result_CXCR4hiCD27midIGMLhi)
barplot(enrich_result_CXCR4hiCD27midIGMLhi)

CXCR4hiMIR155HGhiIGMhi <- ClosestFeature(data, regions = top$CXCR4hiMIR155HGhiIGMhi)
# Perform pathway enrichment
enrich_result_CXCR4hiMIR155HGhiIGMhi <- enrichGO(gene = CXCR4hiMIR155HGhiIGMhi$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.5, qvalueCutoff = 0.5)
# Visualize the enriched pathways
dotplot(enrich_result_CXCR4hiMIR155HGhiIGMhi)
barplot(enrich_result_CXCR4hiMIR155HGhiIGMhi)

CXCR4loCD24hiCLL <- ClosestFeature(data, regions = top$CXCR4loCD24hiCLL)
# Perform pathway enrichment
enrich_result_CXCR4loCD24hiCLL <- enrichGO(gene = CXCR4loCD24hiCLL$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.06, qvalueCutoff = 0.06)
# Visualize the enriched pathways
dotplot(enrich_result_CXCR4loCD24hiCLL)
barplot(enrich_result_CXCR4loCD24hiCLL)

CXCR4loCLL <- ClosestFeature(data, regions = top$CXCR4loCLL)
# Perform pathway enrichment
enrich_result_CXCR4loCLL <- enrichGO(gene = CXCR4loCLL$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.3, qvalueCutoff = 0.3)
# Visualize the enriched pathways
dotplot(enrich_result_CXCR4loCLL)
barplot(enrich_result_CXCR4loCLL)

MZB1hiIGHMhiXBP1hi <- ClosestFeature(data, regions = top$MZB1hiIGHMhiXBP1hi)
# Perform pathway enrichment
enrich_result_MZB1hiIGHMhiXBP1hi <- enrichGO(gene = MZB1hiIGHMhiXBP1hi$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.3, qvalueCutoff = 0.3)
# Visualize the enriched pathways
dotplot(enrich_result_MZB1hiIGHMhiXBP1hi)
barplot(enrich_result_MZB1hiIGHMhiXBP1hi)


## Heatmap w/ DE-genes from ZAP+ and ZAP-
## Get top 100 from ZAP+/ZAP-
data@active.assay = 'RNA'
data <- SetIdent(data, value = as.factor(data$ZAP_expression))

DE_genes_sample <- FindMarkers(object = data, ident.1 = 'ZAP+', ident.2 = 'ZAP-') # Run differential expression analysis
DE_genes_sample = DE_genes_sample[order(-DE_genes_sample[,2], DE_genes_sample[,1] ),]
write.csv(DE_genes_sample, file = "Documents/Documents/after_meeting/data.NOCD3/ATAC/DE_genes_multiome.csv", row.names = T)

DE_genes_sample %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 1.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_genes_sample_pos # aca tengo todos los ZAP70+ genes
write.csv(DE_genes_sample_pos, file = "Documents/Documents/after_meeting/data.NOCD3/ATAC/DE_genes_pos.csv", row.names = T)

DE_genes_sample %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC < -1.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_genes_sample_neg # aca tengo todos los ZAP70- genes
write.csv(DE_genes_sample_neg, file = "Documents/Documents/after_meeting/data.NOCD3/ATAC/DE_genes_neg.csv", row.names = T)

DE_genes_sample_pos = read.csv('Documents/Documents/after_meeting/data.NOCD3/ATAC/DE_genes_pos.csv')
DE_genes_sample_neg = read.csv('Documents/Documents/after_meeting/data.NOCD3/ATAC/DE_genes_neg.csv')

# ZAP70+ 

# Entrez IDs
ZAP_entrez_pos <- bitr(DE_genes_sample_pos$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_sample_pos <- enrichGO(gene = ZAP_entrez_pos$ENTREZID,
                                          OrgDb = org.Hs.eg.db,
                                          ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                          pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                          pvalueCutoff = 0.05,
                                          qvalueCutoff = 0.05)
# Visualization
barplot(enrichment_results_sample_pos)
dotplot(enrichment_results_sample_pos)


# ZAP70-

# Entrez IDs
ZAP_entrez_neg <- bitr(DE_genes_sample_neg$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_sample_neg <- enrichGO(gene = ZAP_entrez_neg$ENTREZID,
                                          OrgDb = org.Hs.eg.db,
                                          ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                          pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                          pvalueCutoff = 0.05,
                                          qvalueCutoff = 0.05)
# Visualization
barplot(enrichment_results_sample_neg)
dotplot(enrichment_results_sample_neg)

####################
# Barplots per groups

proportions <- as.data.frame((table(data$raw_annotation_clusters, data$prognosis) / rowSums(table(data$raw_annotation_clusters, data$prognosis)))*100)
proportions$Freq = sort(proportions$Freq)

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Create proportions data frame
proportions <- as.data.frame((table(data$raw_annotation_clusters, data$prognosis) / rowSums(table(data$raw_annotation_clusters, data$prognosis)))*100)
ggplot(proportions, aes(x = factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Type of prognosis per clusters",
       x = "Clusters",
       y = "Percentage",
       fill = "Prognosis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  scale_fill_manual(values = cols_group)


proportions <- as.data.frame((table(data$orig.ident, data$prognosis) / rowSums(table(data$orig.ident, data$prognosis)))*100)
ggplot(proportions, aes(x = factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Type of prognosis per clusters",
       x = "Clusters",
       y = "Percentage",
       fill = "Prognosis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  scale_fill_manual(values = cols_group)


proportions <- as.data.frame((table(data$ZAP_expression, data$prognosis) / rowSums(table(data$ZAP_expression, data$prognosis)))*100)
ggplot(proportions, aes(x = factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Type of prognosis per clusters",
       x = "Clusters",
       y = "Percentage",
       fill = "Prognosis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  scale_fill_manual(values = cols_group)

library(ggplot2)
library(reshape2)

proportions <- as.data.frame((table(data$ZAP_expression, data$prognosis, data$raw_annotation_clusters)/rowSums(table(data$ZAP_expression, data$prognosis, data$raw_annotation_clusters)))*100)

long_proportions <- melt(proportions, varnames = c("ZAP_expression", "prognosis", "raw_annotation_clusters"), value.name = "proportion")

ggplot(long_proportions, aes(x = Var3, y = proportion, fill = Var1)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Var2) +
  labs(x = "Raw Annotation Clusters", y = "Proportion (%)", fill = "ZAP Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
  scale_fill_manual(values = cols_zap)

## By sample, prognosis, ZAP70 expression and clusters
proportions <- as.data.frame((table(data$orig.ident, data$ZAP_expression, data$prognosis, data$raw_annotation_clusters) / rowSums(table(data$orig.ident, data$ZAP_expression, data$prognosis, data$raw_annotation_clusters)))*100)
long_proportions <- melt(proportions, varnames = c("Sample","ZAP_expression", "prognosis", "raw_annotation_clusters"), value.name = "proportion")

ggplot(long_proportions, aes(x = Var4, y = proportion, fill = interaction(Var3, Var2))) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Var1, scales = "free_x") +
  labs(x = "Cell type raw annotation", y = "Proportion (%)", fill = "Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  scale_fill_manual(values = cols_cond)

### Now that we have all the clusters annotated, we understand the bibliography about ZAP70 and the prognosis of CLL patients related we can generate a signature of genes that are coexpressed with ZAP70 in the bibliography to be able to localize them in our data

features <- c('ZAP70','CD38', 'CD27', 'CXCR4' , 'IGHM', 'MZB1', 'XBP1', 'IGHM','XBP1', 'HDAC9')

agus <- FetchData(data, vars = features, slot = "scale.data")
cluster_info <- data@meta.data$new_annotation

# Aggregate data by clusters
data_aggregated <- aggregate(. ~ cluster_info, agus, mean)
rownames(data_aggregated) <- data_aggregated$cluster_info
data_matrix <- as.matrix(data_aggregated[, -1])  # Remove the cluster column

# Get the order of rows based on cluster labels
cluster_order <- order(data_aggregated$cluster_info, decreasing = T)

# Reorder the rows of data_matrix
data_matrix_ordered <- data_matrix[cluster_order, ]

# Plot the heatmap with ordered rows
pheatmap(data_matrix_ordered,
         cluster_rows = FALSE,  # Since rows are already ordered
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         fontsize_col = 10,
         angle_col = 45, scale = 'column')

features <- unique(features)
DotPlot(data, features = features, group.by = 'new_annotation', cols = cols_sign, dot.scale = 10, col.min = -1.5, col.max = 0.5) + RotatedAxis() 

####################

clusters <- c("CXCR4loIGHMlo", "CXCR4midIGHMlo", "CXCR4loTXNRD1hi","CXCR4loIGHMmidMAST4hi", "CXCR4midIGHMloNOTCH2mid", "CXCR4hiIGHMhi|CD27hiMZB1hiXBP1hi")
data$new_annotation <- factor(data$new_annotation, levels = clusters)

## FINAL PLOTS
cols_anno <- c("CXCR4midIGHMlo"="#287271", "CXCR4loIGHMlo"="#2a9d8f","CXCR4loTXNRD1hi"="#ffcd61","CXCR4hiIGHMhi|CD27hiMZB1hiXBP1hi"="#cf6348","CXCR4midIGHMloNOTCH2mid"="#ee8959", "CXCR4loIGHMmidMAST4hi"="#f7bb8c")
cols_group <- c("good.prognosis"="#a3b18a", "bad.prognosis"="#bc4749")
cols_zap <- c("ZAP-"="#eee3ab", "ZAP+"="#80a1c1")
cols_sign <- c("yellow", "purple")
cols_sam <- c("S1"="#89b0ae", "S2"="#a7c957", "S3"="#5d2a42", "S4"="#ffcb77")
cols_cond <- c('bad.prognosis.ZAP-'='#ece4b7', 'bad.prognosis.ZAP+'='#EEAA49', 'good.prognosis.ZAP-' = '#c1dbb3', 'good.prognosis.ZAP+'='#63AA81')

# basic umaps
DimPlot(data, reduction = 'umap.harmony.no.lymph6', group.by = 'new_annotation', label = F, cols = cols_anno)
DimPlot(data, reduction = 'umap.harmony.no.lymph6', group.by = 'prognosis', label = F, cols = cols_group, order = T)
DimPlot(data, reduction = 'umap.harmony.no.lymph6', group.by = 'ZAP_expression', label = F, cols = cols_zap, order = T)
DimPlot(data, reduction = 'umap.harmony.no.lymph6', group.by = 'orig.ident', label = F, cols = cols_sam)

# signature
FeaturePlot(data1, reduction = "umap.harmony.no.lymph6", features = names(signature), min.cutoff = 0.1, cols = cols_sign, order = T)
# ZAP70
FeaturePlot(data1, reduction = "umap.harmony.no.lymph6", features = names(signature), min.cutoff = 0.1, cols = cols_zap, order = T)


DimPlot(data, reduction = 'umap.harmony.no.lymph6', group.by = 'new_annotation', label = F, cols = cols_anno, split.by = 'orig.ident')

proportions <- as.data.frame((table(my_data$orig.ident, my_data$new_annotation) / rowSums(table(my_data$orig.ident, my_data$new_annotation)))*100)
Clusters_samples <- ggplot(proportions, aes(x = factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Proportion clusters per samples",
       x = "Samples",
       y = "Clusters",
       fill = "Expression") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  scale_fill_manual(values = cols_anno)
