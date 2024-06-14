## Samples reanalysis removing the T cells. (after the first meeting with Marta)

# 1. Sacar las T cells de mi data 
# 2. Hay dos pacientes con mayor expresion (2/3) que los otros dos (1/4) 
# 3. Hacer label transfere (todo y por samples) from the RT data a mi data (matchscore)
# 4. Analizar las leukemic clouds

## My data
data = readRDS('Documents/Documents/MULTIOME/data_oficial.rds') 

# Remove the T cells clusters from my data

t_cell_types <- c("CD4+ T cells", "CD8+ Cytotoxic T cells") # identify CD4+ T cells and CD8+ Cytotoxic T cells
t_cells <- which(data@meta.data$cell_type %in% t_cell_types) # identify T cells (positions)
data_no_Tcells <- subset(data, cells = -t_cells) # subset to remove T cells

data_no_Tcells = readRDS('after_meeting/data_N_Tcells.rds') # read
saveRDS(data_no_Tcells, file = 'after_meeting/data_N_Tcells.rds') # save

table(data_no_Tcells$cell_type)

# Visualize the different metrics with violin plots
Idents(data_no_Tcells) = data_no_Tcells$orig.ident
VlnPlot(data_no_Tcells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

dim(data_no_Tcells)
# 36601 21036

## 1. CLUSTER THE CELLS

data_no_Tcells <- FindNeighbors(data_no_Tcells, dims = 1:20) #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.

data_no_Tcells <- FindClusters(data_no_Tcells, resolution = 0.5) #Partition the KNN graph into clusters.

table(data_no_Tcells@active.ident)

#.  0    1    2    3    4    5    6    7    8 
#. 5814 5161 4834 1820 1324  971  680  244  188 

head(Idents(data_no_Tcells), 5)

## 2. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

data_no_Tcells <- RunUMAP(data_no_Tcells, dims = 1:20, reduction.name = 'umap.rna')


# Samples distribution
DimPlot(data_no_Tcells, group.by = 'sample', reduction = 'umap.rna', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))

# Clusters distribution
DimPlot(data_no_Tcells, group.by = 'seurat_clusters', reduction = 'umap.rna', label = T, cols = c('mediumseagreen','mediumaquamarine','seagreen','cornflowerblue','darkseagreen', 'royalblue','lightblue', '#FFD700', 'salmon1'))
VlnPlot(data_no_Tcells, group.by = 'seurat_clusters', sort = T,features = 'nFeature_RNA', log = T, cols = c('mediumseagreen','mediumaquamarine','seagreen','cornflowerblue','darkseagreen', 'royalblue','lightblue', '#FFD700', 'salmon1'))

# Cell type distribution (no me tengo que quedar con estas anotaciones ya que son muy generales, vamos a tratar de interiorizarnos mas en estos clusters y entender su genetic profile)
DimPlot(data_no_Tcells, group.by = 'cell_type', reduction = 'umap.rna', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen','cornflowerblue', 'pink', 'royalblue', 'lightblue', 'palevioletred','#FFD700'))

## 3. FINDING DIFFERENTIALLY EXPRESSED FEATURES (cluster biomarkers)
Idents(data_no_Tcells) = data_no_Tcells$seurat_clusters
markers <- FindAllMarkers(data_no_Tcells, only.pos = T)

# Save my markers matrix
#write.csv(markers, file = "after_meeting/markers_seurat_clusters.csv", row.names = T)
markers = read.csv('Documents/Documents/MULTIOME/markers_multiome.csv')

markers %>%
  group_by(cluster) %>%
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

DoHeatmap(data_no_Tcells, features = top10$gene) + NoLegend()

library(matchSCore2)

top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=100)
top = as.data.frame(top)
write.csv(top, file = "after_meeting/top100markers_seurat_clusters.csv", row.names = T)
# write.csv(top, file = "Documents/Documents/MULTIOME/top50markers_multiome.csv", row.names = T)
# top = read.csv('Documents/Documents/MULTIOME/top100markers_multiome.csv')


## CLUSTER BY CLUSTER ANALYSIS

# Cluster 0
DimPlot(data_no_Tcells, group.by = "seurat_clusters", cols = c('mediumseagreen','gray','gray','gray','gray', 'gray','gray', 'gray', 'gray'))

head(top$X0, 10)
# "IFNG-AS1" "RPS26" "IGKC" "MAST4" "JCHAIN" "TFEC" "HLA-DRB1" "FUT8" "TMSB4X" "IGHD"  

FeaturePlot(data_no_Tcells, features = head(top$X0, 12), order = T) # location in the clusters

# Cluster 1
DimPlot(data_no_Tcells, group.by = "seurat_clusters", cols = c('gray','mediumaquamarine','gray','gray','gray', 'gray','gray', 'gray', 'gray'))

head(top$X1, 10)
# "IGLC2" "IGLC3" "IGHM" "EIF1" "UBC" "RPL10" "RPS18" "RPS11" "EEF1A1" "CXCR4" 

FeaturePlot(data_no_Tcells, features = head(top$X1, 12), order = T) # location in the clusters

# Cluster 2
DimPlot(data_no_Tcells, group.by = "seurat_clusters", cols = c('gray','gray','seagreen','gray','gray', 'gray','gray', 'gray', 'gray'))

head(top$X2, 10)
# "AC018742.1" "PCDH9" "TSHZ2" "SDK2" "AC068051.1" "TRIO" "CLNK" "ADK" "UFL1-AS1"   "MGAT5"

FeaturePlot(data_no_Tcells, features = head(top$X2, 12), order = T) # location in the clusters

# Cluster 3
DimPlot(data_no_Tcells, group.by = "seurat_clusters", cols = c('gray','gray','gray','cornflowerblue','gray', 'gray','gray', 'gray', 'gray'))

head(top$X3, 10)
# "IGLC2"  "IGLC3"  "RBM20"  "CCDC50" "TRPS1"  "DLG2"   "HSPA5"  "GPM6A"  "GRIK3"  "TESC" 

FeaturePlot(data_no_Tcells, features = head(top$X3, 12), order = T) # location in the clusters

# Cluster 4
DimPlot(data_no_Tcells, group.by = "seurat_clusters", cols = c('gray','gray','gray','gray','darkseagreen', 'gray','gray', 'gray', 'gray'))

head(top$X4, 10)
# "TEX14" "AC253572.2" "AC005580.1" "MT-ND4L" "FOSB" "DMD" "ENPP2" "AC007952.4" "JUND" "JUN" 

FeaturePlot(data_no_Tcells, features = head(top$X4, 12), order = T) # location in the clusters

# Cluster 5
DimPlot(data_no_Tcells, group.by = "seurat_clusters", cols = c('gray','gray','gray','gray','gray', 'royalblue','gray', 'gray', 'gray'))

head(top$X5, 10)
# "IFNG-AS1"  "MAST4" "LINC01374" "JCHAIN" "ALG1L9P" "MVB12B" "SLC44A5" "LIMCH1" "HIF1A" "CAMK4"

FeaturePlot(data_no_Tcells, features = head(top$X5, 12), order = T) # location in the clusters

# Cluster 6
DimPlot(data_no_Tcells, group.by = "seurat_clusters", cols = c('gray','gray','gray','gray','gray', 'gray','lightblue', 'gray', 'gray'))

head(top$X6, 10)
# "AC018742.1" "PCDH9" "UFL1-AS1" "CERS6" "NRCAM" "KLF7" "PTK2" "AL390957.1" "TDRD15" "LINC00640"

FeaturePlot(data_no_Tcells, features = head(top$X6, 12), order = T) # location in the clusters

# Cluster 7
DimPlot(data_no_Tcells, group.by = "seurat_clusters", cols = c('gray','gray','gray','gray','gray', 'gray','gray', '#FFD700', 'gray'))

head(top$X7, 10)
# "DPYD" "AOAH" "PLXDC2" "MYO1F" "NAMPT" "ARHGAP26" "LRMDA" "SLC8A1" "SLCO3A1"  "VCAN" 

FeaturePlot(data_no_Tcells, features = head(top$X7, 12), order = T) # location in the clusters

# Cluster 8
DimPlot(data_no_Tcells, group.by = "seurat_clusters", cols = c('gray','gray','gray','gray','gray', 'gray','gray', 'gray', 'salmon1'))

head(top$X8, 10)
# "PRKCH"   "CD247"   "BCL11B"  "SYTL3"   "TNFAIP3" "STAT4"   "PDE3B"   "CCL5"    "FYB1"    "NKG7" 

FeaturePlot(data_no_Tcells, features = head(top$X8, 12), order = T) # location in the clusters

# Dotplots/heatmap with interesting markers from CLL and the RT paper

specific_markers = c("IFNG-AS1", "JCHAIN", "IGLC2", "IGLC3", "AC018742.1", "UFL1-AS1", "GRIK3", "AC005580.1", "ENPP2", "MVB12B", "TDRD15", "VCAN", "CD14", "IL7R", "GNLY")
DotPlot(data_no_Tcells, features = specific_markers) + RotatedAxis()


## ZAP70 exploration

# Create a column in the metadata for ZAP expression, here we will evaluate if a cell is ZAP positive or not
data_no_Tcells$ZAP_expression <- ifelse(data_no_Tcells$RNA$counts['ZAP70', ] >= 1, "ZAP+", "ZAP-")

table(data_no_Tcells$ZAP_expression)
#  ZAP-  ZAP+ 
# 20151  885

library(ggplot2)
proportions <- as.data.frame((table(data_no_Tcells$seurat_clusters, data_no_Tcells$ZAP_expression) / rowSums(table(data_no_Tcells$seurat_clusters, data_no_Tcells$ZAP_expression)))*100)

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

# EXPRESSION
FeaturePlot(data_no_Tcells, features = 'ZAP70', order  = TRUE, cols = c("grey", "red"), pt.size = 0.5, reduction = 'umap.rna')

VlnPlot(data_no_Tcells, features = 'ZAP70', group.by = 'sample', sort = T, log = T)

# COLOR CELLS THAT ARE ZAP+ 
DimPlot(data_no_Tcells, group.by = "ZAP_expression", cols = c("lightgray", "red"), pt.size = 0.1, order = T, reduction = 'umap.rna')

## RELATED MARKERS w/ ZAP70 from the RT dataset

FeaturePlot(data_no_Tcells, features = c('CXCR4', 'CD27', 'TP53INP1'), order = T, max.cutoff = 'q90', label = T) # location in the clusters

p1 = FeaturePlot(data_no_Tcells, features = 'CXCR4', order = T, max.cutoff = 'q90', label = T) # location in the clusters
p2 = VlnPlot(data_no_Tcells, features = 'CXCR4', group.by = 'seurat_clusters', sort = T, cols = c('salmon1', 'cornflowerblue', 'mediumaquamarine','lightblue', 'seagreen', 'darkseagreen','royalblue','mediumseagreen', '#FFD700' ))
p1*p2

p1 = FeaturePlot(data_no_Tcells, features = 'CD27', order = T, max.cutoff = 'q90', label = T) # location in the clusters
p2 = VlnPlot(data_no_Tcells, features = 'CD27', group.by = 'seurat_clusters', sort = T, cols = c('salmon1', 'royalblue', 'mediumseagreen','cornflowerblue', 'mediumaquamarine', 'darkseagreen','lightblue','seagreen', '#FFD700' ))
p1*p2

p1 = FeaturePlot(data_no_Tcells, features = 'TP53INP1', order = T, max.cutoff = 'q90', label = T) # location in the clusters
p2 = VlnPlot(data_no_Tcells, features = 'TP53INP1', group.by = 'seurat_clusters', sort = T, cols = c('royalblue', 'lightblue', 'mediumseagreen','salmon1', 'seagreen', '#FFD700','darkseagreen','cornflowerblue', 'mediumaquamarine' ))
p1*p2

# Coexpression with ZAP70
FeaturePlot(data_no_Tcells, features = c("ZAP70", "CXCR4"), blend = TRUE, order = T)
FeaturePlot(data_no_Tcells, features = c("ZAP70", "CD27"), blend = TRUE, order = T)
FeaturePlot(data_no_Tcells, features = c("ZAP70", "TP53INP1"), blend = TRUE, order = T)

VlnPlot(data_no_Tcells, features = c('CXCR4', 'ZAP70', 'CD27', 'TP53INP1'), group.by = 'cell_type', sort = T)

VlnPlot(data_no_Tcells, features = c('CXCR4', 'CD27'), group.by = 'ZAP_expression')

VlnPlot(data_no_Tcells, features = c('CXCR4', 'CD27', 'MIR155HG', 'TP53INP1', 'CCND2'), group.by = 'seurat_clusters', split.by = 'ZAP_expression')
# rosa ZAP70-
# azul ZAP70+

## Chronic Lymphocytic leukemia markers (different levels of expression)
VlnPlot(data_no_Tcells, features = c('CXCR4', 'CD27', 'MIR155HG'), group.by = 'seurat_clusters')

## Richter transformation markers (different levels of expression)
VlnPlot(data_no_Tcells, features = c('CCND2', 'MIR155HG', 'TP53INP1'), group.by = 'seurat_clusters')

## RT intraclonal heterogeneity highly expressed markers
VlnPlot(data_no_Tcells, features = c('MKI67', 'PCNA'), group.by = 'seurat_clusters')

## IDEAS
FeaturePlot(data_no_Tcells, features = c("CXCR4", "CD27", "MIR155HG", "ZAP70"), order = T)

VlnPlot(data_no_Tcells, features = c("CXCR4", "CD27", "MIR155HG", "ZAP70"), group.by = 'seurat_clusters', sort = T)

VlnPlot(data_no_Tcells, features = c("MIR155HG", "ZAP70"), group.by = 'seurat_clusters', sort = T)

FeaturePlot(data_no_Tcells, features = c("CXCR4", "CD27"), blend = TRUE, order = T)

VlnPlot(data_no_Tcells, features = c("MIR155HG", "ZAP70"), group.by = 'seurat_clusters', sort = T)

FeaturePlot(data_no_Tcells, features = c("ZAP70", "IGVH"), blend = TRUE, order = T)

## External dataset
data_rt = readRDS('RT/data_rt.rds')
'''
### Training of the model  

Idents(data_rt) = data_rt$annotation_final
markers <- FindAllMarkers(data_rt, only.pos = T)
#markers = read.csv('RT/top100markers_RT.csv', row.names = F)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=100)

scaled_ref <- ScaleData(data_rt, features = unlist(top))

## Generate the model 
model <- train_model(scale.data = scaled_ref,clus = data_rt$annotation_final, gene_cl.ref = top,prop = 0.75)

## Cell projection in my data
out <- identity_map(scale.data = data_no_Tcells@assays$RNA$scale.data,model = model,top)

### cell identities
ids <- out$ids 

# Save the model
saveRDS(model, file = "after_meeting/model_RT_annotations.rds")
model = readRDS(file = 'after_meeting/model_RT_annotations.rds')
'''
data_no_Tcells$matchScore_RT_annotation <- ids # add to my metadata the annotation of my model trained with the RT dataset

table(data_no_Tcells$matchScore_RT_annotation)

# CCND2hi RT = 6
# CCND2lo RT = 7 

# CLL = 4

# CXCR4hiCD27lo = 10938
# CXCR4loCD27hi = 2796
# CXCR4loCD27hi RT = 2

# MIR155HGhi RT = 1137

# MZB1hiIGHMhiXBP1hi = 1
# RT proliferative = 5214
# TP53INP1hi RT = 366

# unclassified = 565

p1 = DimPlot(data_no_Tcells, group.by = 'matchScore_RT_annotation', reduction = 'umap.rna', label = F, cols = c('gray', 'gray', 'gray', 'yellow', 'red', 'gray', 'green', 'gray', 'pink', 'orange', 'blue'), order = T)
p2 = DimPlot(data_no_Tcells, group.by = 'seurat_clusters', reduction = 'umap.rna', order = T, label = T, cols = c('mediumseagreen','mediumaquamarine','seagreen','cornflowerblue','darkseagreen', 'royalblue','lightblue', '#FFD700', 'salmon1'))
p1 + p2

table(data_no_Tcells$matchScore_RT_annotation, data_no_Tcells$seurat_clusters)
# important annotations: CXCR4hiCD27lo, CXCR4loCD27hi, MIR155HGhi RT, RT proliferative, TP53INP1hi RT, unclassified

table(data_no_Tcells$matchScore_RT_annotation == 'CXCR4hiCD27lo', data_no_Tcells$seurat_clusters)
#    0    1    2    3    4    5    6    7    8
# 4665  603 3527  105  910  662  445   20    1

table(data_no_Tcells$matchScore_RT_annotation == 'CXCR4loCD27hi', data_no_Tcells$seurat_clusters)
#    0    1    2    3    4    5    6    7    8
#  754  474  754  209   76  281  196   25   27

table(data_no_Tcells$matchScore_RT_annotation == 'MIR155HGhi RT', data_no_Tcells$seurat_clusters)
#    0    1    2    3    4    5    6    7    8
#   108  755   65   90   41    1    4   55   18

table(data_no_Tcells$matchScore_RT_annotation == 'RT proliferative', data_no_Tcells$seurat_clusters)
#    0    1    2    3    4    5    6    7    8
#   230 2894  436 1125  235   21   28  140  105

table(data_no_Tcells$matchScore_RT_annotation == 'TP53INP1hi RT', data_no_Tcells$seurat_clusters)
#    0    1    2    3    4    5    6    7    8
#   14  128    9  158   34    3    2    0   18

table(data_no_Tcells$matchScore_RT_annotation == 'unclassified', data_no_Tcells$seurat_clusters)
#    0    1    2    3    4    5    6    7    8
#   42  300   40  129   25    2    5    4   18


## barplot to observe the proportion of RT_annotations per cluster per sample

table(data_no_Tcells$matchScore_RT_annotation) # RT annotations
table(data_no_Tcells$sample) # samples
table(data_no_Tcells$seurat_clusters) # clusters

proportions = as.data.frame((table(data_no_Tcells$sample, data_no_Tcells$seurat_clusters, data_no_Tcells$matchScore_RT_annotation) / rowSums(table(data_no_Tcells$sample, data_no_Tcells$seurat_clusters, data_no_Tcells$matchScore_RT_annotation)))*100) 


ggplot(proportions, aes(x = as.factor(Var2), y = Freq, fill = Var3)) +
  geom_bar(stat = "identity", position = "fill") +  # Change position to "fill"
  labs(x = "Clusters", y = "Proportion", fill = "matchScore_RT_annotation") +
  theme_bw() +
  theme(legend.position = "bottom") +   
  scale_fill_manual(values = c("cornflowerblue", "lightblue","gold","darkseagreen", "palegreen4", "red4","palevioletred","orange","palegoldenrod","pink", "gray"))

#+ facet_grid(~Var1)

########

Idents(data_no_Tcells) = data_no_Tcells$matchScore_RT_annotation
markers <- FindAllMarkers(data_no_Tcells, only.pos = T)
#write.csv(markers, file = 'after_meeting/markers_RT_annotation.csv')

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

DoHeatmap(data_no_Tcells, features = top10$gene) + NoLegend()

top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=100)
top = as.data.frame(top)
#write.csv(top, file = "after_meeting/top100markers_RT_annotation.csv", row.names = T)


VlnPlot(data_no_Tcells, features = 'ZAP70', sort = T, group.by = 'matchScore_RT_annotation', log = T)

Idents(data_no_Tcells) = as.factor(data_no_Tcells$matchScore_RT_annotation)

DotPlot(data_no_Tcells, features = c('CXCR4', 'CD27', 'CCND2', 'MIR155HG', 'TP53INFP1', 'ZAP70'), col.min = 0) + RotatedAxis()

table(data_no_Tcells$matchScore_RT_annotation)


table(data_no_Tcells$matchScore_RT_annotation, data_no_Tcells$seurat_clusters)
VlnPlot(data_no_Tcells, features = 'ZAP70', sort = T, group.by = 'seurat_clusters', split.by = 'matchScore_RT_annotation',log = T)


VlnPlot(data_no_Tcells, features = c('CD38', 'CXCR4', 'CD27', 'STAT3', 'CD5'), group.by = 'matchScore_RT_annotation', sort = T, log = T)

DotPlot(data_no_Tcells, features = c('CD38', 'CXCR4', 'CD27', 'STAT3', 'CD5'), group.by = 'ZAP_expression')

# Per sample (my dataset)
# IMPORTANT!! the top variable is from the RT dataset

# Here we will applied the previouly generated model in order to infere the celltypes of the RT dataset in my own data (sample by sample) 
# We want to observe how the labels from the RT dataset are distributed across our samples individually

S1 = readRDS('Documents/Documents/S1.rds')
out_S1 <- identity_map(scale.data = S1@assays$RNA$scale.data,model = model,top)
ids_S1 <- out_S1$ids 

table(ids_S1)

# CCND2lo RT      CLL    CXCR4hiCD27lo    CXCR4loCD27hi 
#       2          5         2943             1368 

# MIR155HGhi RT    RT proliferative    TP53INP1hi RT     unclassified 
#      664              852               25                 108 

S1$matchScore_RT_annotation <- ids_S1 # add to my metadata the annotation of my model trained with the RT dataset

p1 = DimPlot(S1, group.by = 'matchScore_RT_annotation', reduction = 'umap', label = F, cols = c('orange','gold','cornflowerblue','palegreen4', 'red4','palevioletred', 'palegoldenrod', 'gray'), order = T)
p2 = DimPlot(S1, reduction = "umap", label = T, pt.size = 0.5,  cols = c('#F4D35E','violet','pink','palevioletred','mediumseagreen','mediumaquamarine','darkseagreen','seagreen')) + NoLegend()
p1 + p2

saveRDS(S1, file = 'S1.rds')

############

S2 = readRDS('S2.rds')

out_S2 <- identity_map(scale.data = S2@assays$RNA$scale.data,model = model,top)
ids_S2 <- out_S2$ids 

table(ids_S2)
#   CCND2hi RT       CCND2lo RT     CLL    CXCR4hiCD27lo    CXCR4loCD27hi 
#        2                1          5         4202             1318 

# CXCR4loCD27hi RT    MIR155HGhi RT   RT proliferative    TP53INP1hi RT     unclassified 
#           1              609             1480               54              148 

S2$matchScore_RT_annotation <- ids_S2 # add to my metadata the annotation of my model trained with the RT dataset

p1 = DimPlot(S2, group.by = 'matchScore_RT_annotation', reduction = 'umap', label = F, cols = c('orange','gold','mediumaquamarine','palegreen4', 'cornflowerblue','red4', 'palevioletred', 'palegoldenrod','pink', 'gray'), order = T)
p2 = DimPlot(S2, reduction = "umap", label = T, pt.size = 0.5, cols = c('violet','pink','mediumseagreen','mediumaquamarine','darkseagreen','seagreen', 'green4', '#F4D35E','palevioletred3', 'palevioletred1')) + NoLegend()
p1 + p2

saveRDS(S2, file = 'S2.rds')

############

S3 = readRDS('S3.rds')
out_S3 <- identity_map(scale.data = S3@assays$RNA$scale.data,model = model,top)
ids_S3 <- out_S3$ids 

table(ids_S3)

S3$matchScore_RT_annotation <- ids_S3 # add to my metadata the annotation of my model trained with the RT dataset

p1 = DimPlot(S3, group.by = 'matchScore_RT_annotation', reduction = 'umap', label = F, cols = c('orange','gold','cornflowerblue','palegreen4', 'mediumaquamarine','palevioletred', 'red4', 'gray'), order = T)
p2 = DimPlot(S3, reduction = "umap", label = T, pt.size = 0.5, cols = c('violet','pink','mediumseagreen','mediumaquamarine','darkseagreen','seagreen', 'green4', '#F4D35E')) + NoLegend()

p1 + p2

saveRDS(S3, file = 'S3.rds')

############

S4 = readRDS('S4.rds')
out_S4 <- identity_map(scale.data = S4@assays$RNA$scale.data,model = model,top)
ids_S4 <- out_S4$ids 

table(ids_S4)

S4$matchScore_RT_annotation <- ids_S4 # add to my metadata the annotation of my model trained with the RT dataset

p1 = DimPlot(S4, group.by = 'matchScore_RT_annotation', reduction = 'umap', label = F, cols = c('orange','gold','mediumaquamarine','palegreen4', 'cornflowerblue','red4', 'palevioletred', 'palegoldenrod','pink', 'gray'), order = T)
p2 = DimPlot(S4, reduction = "umap", label = T, pt.size = 0.5, cols = c('violet','pink','mediumseagreen','mediumaquamarine','darkseagreen','seagreen')) + NoLegend()

p1 + p2

saveRDS(S4, file = 'S4.rds')

#######################
# SAMPLE 1 - reanalysis

# remove the Tcells and myeloid cluster from the sample
S1_t_cells <- c("Cytotoxic NK","CD4+ T cells", "CD8+ Cytotoxic T cells","CD14+ Mono") 
t_cells <- which(S1@meta.data$cell_type %in% S1_t_cells) # identify T cells (positions)
S1_no_Tcells <- subset(S1, cells = -t_cells) # subset to remove T cells
saveRDS(S1_no_Tcells, file = 'after_meeting/S1_no_Tcells.rds')
S1_no_Tcells = readRDS('after_meeting/S1_no_Tcells.rds')

# Visualize the different metrics with violin plots
Idents(S1_no_Tcells) = S1_no_Tcells$orig.ident
VlnPlot(S1_no_Tcells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

dim(S1_no_Tcells)
# 36601 5490

## 1. CLUSTER THE CELLS

S1_no_Tcells <- FindNeighbors(S1_no_Tcells, dims = 1:20) #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.
S1_no_Tcells <- FindClusters(S1_no_Tcells, resolution = 0.5) #Partition the KNN graph into clusters.
table(S1_no_Tcells@active.ident)
#.    0   1    2     3       
#.  2014 1395 1270  811 

head(Idents(S1_no_Tcells), 5)

## 2. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

S1_no_Tcells <- RunUMAP(S1_no_Tcells, dims = 1:20, reduction.name = 'umap.rna')

# Clusters distribution
DimPlot(S1_no_Tcells, group.by = 'seurat_clusters', reduction = 'umap.rna', label = T, cols = c('palegreen4','lightcoral','lightgoldenrod2','lightskyblue3'))
VlnPlot(S1_no_Tcells, group.by = 'seurat_clusters', sort = T,features = 'nFeature_RNA', log = T, cols = c('lightskyblue3','palegreen4','lightgoldenrod2','lightcoral'))

# Cell type distribution (no me tengo que quedar con estas anotaciones ya que son muy generales, vamos a tratar de interiorizarnos mas en estos clusters y entender su genetic profile)
# DimPlot(S1_no_Tcells, group.by = 'cell_type', reduction = 'umap.rna', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen','cornflowerblue', 'pink', 'royalblue', 'lightblue', 'palevioletred','#FFD700'))

table(S1_no_Tcells$matchScore_RT_annotation)

# CCND2lo RT     CLL    CXCR4hiCD27lo    CXCR4loCD27hi    MIR155HGhi RT  RT proliferative    TP53INP1hi RT     unclassified 
#    2            5        2889             1274              503              695               24                 98 

p1 = DimPlot(S1_no_Tcells, group.by = 'matchScore_RT_annotation', reduction = 'umap.rna', label = F, cols = c('gray','gray','seagreen','cornflowerblue','pink','palevioletred', '#FFD700','gray'), order = T)
p2 = DimPlot(S1_no_Tcells, group.by = 'seurat_clusters', reduction = 'umap.rna', order = T, label = T, cols = c('palegreen4','lightcoral','lightgoldenrod2','lightskyblue3'))
p1 + p2

p2 + p1


## barplot to observe the proportion of RT_annotations per cluster 
proportions = as.data.frame((table(S1_no_Tcells$seurat_clusters, S1_no_Tcells$matchScore_RT_annotation) / rowSums(table(S1_no_Tcells$seurat_clusters, S1_no_Tcells$matchScore_RT_annotation)))*100) 

ggplot(proportions, aes(x = as.factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +  # Change position to "fill"
  labs(x = "Clusters sample 1", y = "Proportion", fill = "matchScore_RT_annotation") +
  theme_bw() +
  theme(legend.position = "bottom") +   
  scale_fill_manual(values = c("cornflowerblue", "lightblue","darkseagreen","palegreen4", "palegoldenrod", "orange","palevioletred","gray"))

# ZAP70 EXPRESSION 
FeaturePlot(S1_no_Tcells, features = 'ZAP70', order  = TRUE, cols = c("grey", "red"), pt.size = 0.5, reduction = 'umap.rna')
# by RT annotation
VlnPlot(S1_no_Tcells[, S1_no_Tcells$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'matchScore_RT_annotation', log = TRUE, sort = T, cols = c('cornflowerblue','seagreen'))

#######################
# SAMPLE 2 - reanalysis

# remove the Tcells and myeloid cells from the sample
S2_t_cells <- c("CD4+ T cells", "CD8+ cytotoxic T cells", "T cells (FCGR3A+)", "T cells (GNLY+)", "Myeloid cells")
t_cells <- which(S2@meta.data$cell_type %in% S2_t_cells) # identify T cells (positions)
S2_no_Tcells <- subset(S2, cells = -t_cells) # subset to remove T cells
saveRDS(S2_no_Tcells, file = 'after_meeting/S2_no_Tcells.rds')

# Visualize the different metrics with violin plots
Idents(S2_no_Tcells) = S2_no_Tcells$orig.ident
VlnPlot(S2_no_Tcells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

dim(S2_no_Tcells)
# 36601 7339

## 1. CLUSTER THE CELLS

S2_no_Tcells <- FindNeighbors(S2_no_Tcells, dims = 1:20) #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.
S2_no_Tcells <- FindClusters(S2_no_Tcells, resolution = 0.5) #Partition the KNN graph into clusters.
table(S2_no_Tcells@active.ident)
#.    0    1   2     3    4   
#.  2199 2038 1495 1284  323 

head(Idents(S2_no_Tcells), 5)

## 2. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

S2_no_Tcells <- RunUMAP(S2_no_Tcells, dims = 1:20, reduction.name = 'umap.rna')

# Clusters distribution
DimPlot(S2_no_Tcells, group.by = 'seurat_clusters', reduction = 'umap.rna', label = T, cols = c('palegreen4','lightcoral','lightgoldenrod2','lightskyblue3', 'rosybrown2'))
VlnPlot(S2_no_Tcells, group.by = 'seurat_clusters', sort = T,features = 'nFeature_RNA', log = T, cols = c('rosybrown2','lightgoldenrod2','lightskyblue3','lightcoral', 'palegreen4'))

# Cell type distribution (no me tengo que quedar con estas anotaciones ya que son muy generales, vamos a tratar de interiorizarnos mas en estos clusters y entender su genetic profile)
#DimPlot(S2_no_Tcells, group.by = 'cell_type', reduction = 'umap.rna', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen','cornflowerblue', 'pink', 'royalblue', 'lightblue', 'palevioletred','#FFD700'))

table(S2_no_Tcells$matchScore_RT_annotation)

# CCND2hi RT       CCND2lo RT              CLL    CXCR4hiCD27lo    CXCR4loCD27hi CXCR4loCD27hi RT  MIR155HGhi RT RT proliferative    TP53INP1hi RT     unclassified 
#.    2                1                5             4175             1217                1           352             1415               47              124 

p1 = DimPlot(S2_no_Tcells, group.by = 'matchScore_RT_annotation', reduction = 'umap.rna', label = F, cols = c('gray','gray','gray','seagreen','cornflowerblue','gray', 'pink','palevioletred', '#FFD700', 'gray'), order = T)
p2 = DimPlot(S2_no_Tcells, group.by = 'seurat_clusters', reduction = 'umap.rna', order = T, label = T, cols = c('palegreen4','lightcoral','lightgoldenrod2','lightskyblue3', 'rosybrown2'))
p1 + p2
p2 + p1

## barplot to observe the proportion of RT_annotations per cluster
proportions = as.data.frame((table(S2_no_Tcells$seurat_clusters, S2_no_Tcells$matchScore_RT_annotation) / rowSums(table(S2_no_Tcells$seurat_clusters, S2_no_Tcells$matchScore_RT_annotation)))*100) 

ggplot(proportions, aes(x = as.factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +  # Change position to "fill"
  labs(x = "Clusters sample 2", y = "Proportion", fill = "matchScore_RT_annotation") +
  theme_bw() +
  theme(legend.position = "bottom") +   
  scale_fill_manual(values = c("cornflowerblue", "lightblue","royalblue","darkseagreen","palegreen4", "mediumseagreen", "orange","palegoldenrod","palevioletred","gray"))

# ZAP70 EXPRESSION 
FeaturePlot(S2_no_Tcells, features = 'ZAP70', order  = TRUE, cols = c("grey", "red"), pt.size = 0.5, reduction = 'umap.rna')
# by RT annotation
VlnPlot(S2_no_Tcells[, S2_no_Tcells$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'matchScore_RT_annotation', log = TRUE, sort = T, cols = c('pink','seagreen','#FFD700', 'palevioletred', 'gray', 'cornflowerblue'))

#######################
# SAMPLE 3 - reanalysis

# remove the Tcells and myeloid cluster from the sample
S3_t_cells <- c("CD4+ T cells", "CD8+ Cytotoxic T cells", "Myeloid cells")
t_cells <- which(S3@meta.data$cell_type %in% S3_t_cells) # identify T cells (positions)
S3_no_Tcells <- subset(S3, cells = -t_cells) # subset to remove T cells
saveRDS(S3_no_Tcells, file = 'after_meeting/S3_no_Tcells.rds')

# Visualize the different metrics with violin plots
Idents(S3_no_Tcells) = S3_no_Tcells$orig.ident
VlnPlot(S3_no_Tcells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

dim(S3_no_Tcells)
# 36601 1679

## 1. CLUSTER THE CELLS

S3_no_Tcells <- FindNeighbors(S3_no_Tcells, dims = 1:20) #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.
S3_no_Tcells <- FindClusters(S3_no_Tcells, resolution = 0.5) #Partition the KNN graph into clusters.
table(S3_no_Tcells@active.ident)
#.  0    1   2    3    4     
#. 539  463 305  277  95 

head(Idents(S3_no_Tcells), 5)

## 2. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

S3_no_Tcells <- RunUMAP(S3_no_Tcells, dims = 1:20, reduction.name = 'umap.rna')

# Clusters distribution
DimPlot(S3_no_Tcells, group.by = 'seurat_clusters', reduction = 'umap.rna', label = T, cols = c('palegreen4','lightcoral','lightgoldenrod2','lightskyblue3', 'rosybrown2'))
VlnPlot(S3_no_Tcells, group.by = 'seurat_clusters', sort = T,features = 'nFeature_RNA', log = T, cols = c('rosybrown2','lightgoldenrod2','palegreen4','lightskyblue3', 'lightcoral'))

# Cell type distribution (no me tengo que quedar con estas anotaciones ya que son muy generales, vamos a tratar de interiorizarnos mas en estos clusters y entender su genetic profile)
# DimPlot(S3_no_Tcells, group.by = 'cell_type', reduction = 'umap.rna', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen','cornflowerblue', 'pink', 'royalblue', 'lightblue', 'palevioletred','#FFD700'))

table(S3_no_Tcells$matchScore_RT_annotation)

# CLL    CXCR4hiCD27lo    CXCR4loCD27hi    MIR155HGhi RT RT proliferative    TP53INP1hi RT  unclassified
#  1        962              288               88              277               19             44 

p1 = DimPlot(S3_no_Tcells, group.by = 'matchScore_RT_annotation', reduction = 'umap.rna', label = F, cols = c('gray','seagreen','cornflowerblue','pink','palevioletred','#FFD700', 'gray'), order = T)
p2 = DimPlot(S3_no_Tcells, group.by = 'seurat_clusters', reduction = 'umap.rna', order = T, label = T, cols = c('palegreen4','lightcoral','lightgoldenrod2','lightskyblue3', 'rosybrown2'))
p1 + p2
p2 + p1

## barplot to observe the proportion of RT_annotations per cluster
proportions = as.data.frame((table(S3_no_Tcells$seurat_clusters, S3_no_Tcells$matchScore_RT_annotation) / rowSums(table(S3_no_Tcells$seurat_clusters, S3_no_Tcells$matchScore_RT_annotation)))*100) 

ggplot(proportions, aes(x = as.factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +  # Change position to "fill"
  labs(x = "Clusters", y = "Proportion", fill = "matchScore_RT_annotation") +
  theme_bw() +
  theme(legend.position = "bottom") +   
  scale_fill_manual(values = c("lightblue", "darkseagreen","palegreen4", "palegoldenrod", "orange","palevioletred", "gray"))

# ZAP70 EXPRESSION 
FeaturePlot(S3_no_Tcells, features = 'ZAP70', order  = TRUE, cols = c("grey", "red"), pt.size = 0.5, reduction = 'umap.rna')
# by RT annotation
VlnPlot(S3_no_Tcells[, S3_no_Tcells$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'matchScore_RT_annotation', log = TRUE, sort = T, cols = c('pink','#FFD700','gray', 'palevioletred', 'seagreen', 'cornflowerblue','gray'))

#######################
# SAMPLE 4 - reanalysis

# remove the Tcells from the sample (we couldnt recognise any myeloid cluster)
S4_t_cells <- c("CD4+ T cells", "CD8+ Cytotoxic T cells")
t_cells <- which(S4@meta.data$cell_type %in% S4_t_cells) # identify T cells (positions)
S4_no_Tcells <- subset(S4, cells = -t_cells) # subset to remove T cells
saveRDS(S4_no_Tcells, file = 'after_meeting/S4_no_Tcells.rds')

# Visualize the different metrics with violin plots
Idents(S4_no_Tcells) = S4_no_Tcells$orig.ident
VlnPlot(S4_no_Tcells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

dim(S4_no_Tcells)
# 36601 6199

## 1. CLUSTER THE CELLS

S4_no_Tcells <- FindNeighbors(S4_no_Tcells, dims = 1:20) #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.
S4_no_Tcells <- FindClusters(S4_no_Tcells, resolution = 0.5) #Partition the KNN graph into clusters.
table(S4_no_Tcells@active.ident)
#.  0    1     2    3        
#. 3348 2109  380  362   

head(Idents(S4_no_Tcells), 5)

## 2. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

S4_no_Tcells <- RunUMAP(S4_no_Tcells, dims = 1:20, reduction.name = 'umap.rna')

# Clusters distribution
DimPlot(S4_no_Tcells, group.by = 'seurat_clusters', reduction = 'umap.rna', label = T, cols = c('palegreen4','lightcoral','lightgoldenrod2','lightskyblue3'))
VlnPlot(S4_no_Tcells, group.by = 'seurat_clusters', sort = T,features = 'nFeature_RNA', log = T, cols = c('lightgoldenrod2','lightcoral','lightskyblue3','palegreen4'))

# Cell type distribution (no me tengo que quedar con estas anotaciones ya que son muy generales, vamos a tratar de interiorizarnos mas en estos clusters y entender su genetic profile)
# DimPlot(S4_no_Tcells, group.by = 'cell_type', reduction = 'umap.rna', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen','cornflowerblue', 'pink', 'royalblue', 'lightblue', 'palevioletred','#FFD700'))

table(S4_no_Tcells$matchScore_RT_annotation)

# CCND2hi RT       CCND2lo RT        CLL    CXCR4hiCD27lo    CXCR4loCD27hi    MIR155HGhi RT 
#.  8                4                6             3002             1086              697 

# RT proliferative    TP53INP1hi RT     unclassified 
#.   1237                7              152 

p1 = DimPlot(S4_no_Tcells, group.by = 'matchScore_RT_annotation', reduction = 'umap.rna', label = F, cols = c('gray','gray','gray','seagreen','cornflowerblue','pink', 'palevioletred', '#FFD700', 'gray'), order = T)
p2 = DimPlot(S4_no_Tcells, group.by = 'seurat_clusters', reduction = 'umap.rna', order = T, label = T, cols = c('palegreen4','lightcoral','lightgoldenrod2','lightskyblue3'))
p1 + p2
p2 + p1

## barplot to observe the proportion of RT_annotations per cluster
proportions = as.data.frame((table(S4_no_Tcells$seurat_clusters, S4_no_Tcells$matchScore_RT_annotation) / rowSums(table(S4_no_Tcells$seurat_clusters, S4_no_Tcells$matchScore_RT_annotation)))*100) 

ggplot(proportions, aes(x = as.factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +  # Change position to "fill"
  labs(x = "Clusters", y = "Proportion", fill = "matchScore_RT_annotation") +
  theme_bw() +
  theme(legend.position = "bottom") +   
  scale_fill_manual(values = c("royalblue","cornflowerblue", "lightblue","darkseagreen","palegreen4", "palegoldenrod", "orange","palevioletred", "gray"))

# ZAP70 EXPRESSION 
FeaturePlot(S4_no_Tcells, features = 'ZAP70', order  = TRUE, cols = c("grey", "red"), pt.size = 0.5, reduction = 'umap.rna')
# by RT annotation
VlnPlot(S4_no_Tcells[, S4_no_Tcells$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'matchScore_RT_annotation', log = TRUE, sort = T, cols = c('pink','palevioletred','cornflowerblue','seagreen', 'gray'))

#######################

#ZAP70 proportions per sample

library(ggplot2)

proportions <- as.data.frame((table(S4_no_Tcells$seurat_clusters, S4_no_Tcells$ZAP_expression) / rowSums(table(S4_no_Tcells$seurat_clusters, S4_no_Tcells$ZAP_expression)))*100)
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

#######################
# TO MERGE THE RNA ASSAYS FROM THE FOUR SAMPLES
S1_no_Tcells$orig.ident = 'S1'
S2_no_Tcells$orig.ident = 'S2'
S3_no_Tcells$orig.ident = 'S3'
S4_no_Tcells$orig.ident = 'S4'

# Get the RNA assays of all samples
S1_RNA_assay = DietSeurat(object = S1_no_Tcells, assays = "RNA")
S2_RNA_assay = DietSeurat(object = S2_no_Tcells, assays = "RNA")
S3_RNA_assay = DietSeurat(object = S3_no_Tcells, assays = "RNA")
S4_RNA_assay = DietSeurat(object = S4_no_Tcells, assays = "RNA")

# Merge the four samples
data <- merge(x = S1_RNA_assay,
                  y = list(S2_RNA_assay, S3_RNA_assay, S4_RNA_assay),
                  add.cell.ids = c("S1", "S2", "S3", "S4"))

## Pipeline for RNA data

data@active.assay = 'RNA'
dim(data)
# 36601 20707

head(rownames(data), 5) # genes
head(colnames(data), 5) # cells 

data$object <- 'SeuratObject'

Idents(data) = data$object
# Visualize the different metrics with violin plots
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

hist(data$percent.mt, main = 'Distribution of mitochondrial genes in all samples')
hist(data$nFeature_RNA, main = 'Distribution of the nFeatures in all samples')
plot(density(data$nFeature_RNA)) + abline(v = 250)
plot(density(data$percent.mt)) 

summary(data[['nFeature_RNA']])
summary(data[['percent.mt']])

# Visualize feature-feature relationship with scatterplots
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
# We expect this plot to not have any type of correlation as plotting the depth of the samples with the percentage of mitochondrial RNA. In fact we wan more RNA without mitochondrial. The correlation coefficient is 0
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# We expect this plot to have a strong correlation as we are relating the number of RNA counts (depth) and the number of features meaning the genes. 
# This can be explained as the higher the depth the better the coverage and sensitivity in detecting gene expression making the analysis much more 
# comprehensive analysis, moreover higher numbers of features suggest greater biological complexity and diversity in terms of the genes being expressed.

# Deeper sequencing will result in the detection of a greater number of expressed genes, reflecting the biological complexity of the cells being analyzed.
plot1+plot2

### SUBSET IF NECESSARY ###

summary(data$nFeature_RNA) # It is really interesting to look at the summary of the number of features as it helps us to choose a threshold to remove problematic cells
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#. 504    1072    1287    1542    1617    9806 

summary(data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   2.806   4.012   4.653   5.877  19.966  

### SELECT THOSE CELLS THAT HAVE A GOOD SCORE REGARDING NUMBER OF FEATURES AND THE PERCENTAGE OF MITOCHONDRIAL CONTAMINATION!!!

#data <- subset(data, subset = nFeature_RNA > 800 & percent.mt < 20) # We are going to select the cells we want to be analysed based on the previous results mainly the violin plots from above
#dim(data)

### SAVE ###
# Visualize the different metrics with violin plots
# VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

# Visualize feature-feature relationship with scatterplots
# plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1+plot2

# summary(data@meta.data$nFeature_RNA) # It is really interesting to look at the summary of the number of features as it helps us to choose a threshold to remove problematic cells
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   809    1243    1412    1611    1684    8918 

# summary(data@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2697  4.3312  5.5556  5.7388  7.0505  9.9877 

#plot(density(log10(data@meta.data$nFeature_RNA)))

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

data <- ScaleData(data) # This function scales and centers the expression values for each gene across all cells. 

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
data <- FindClusters(data, resolution = 0.5) #Partition the KNN graph into clusters.
table(data@active.ident)

#.   0    1    2    3    4    5    6    7    8     
#. 5515 5051 5044 1504 1467 1158  417  340  211 


## 8. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

data <- RunUMAP(data, dims = 1:20, reduction.name = 'umap.rna_0.5')

DimPlot(data, group.by = 'orig.ident', reduction = 'umap.rna_0.5', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))
VlnPlot(data[, data$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'orig.ident', log = TRUE, sort = T, cols = c('mediumseagreen','darkseagreen','seagreen','mediumaquamarine'))

p1 = DimPlot(data, group.by = 'matchScore_RT_annotation', reduction = 'umap.rna_0.5', label = F, cols = c('gray','gray','gray','darkseagreen','seagreen', 'gray', 'red', 'yellow', 'black','gray'))
VlnPlot(data[, data$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'matchScore_RT_annotation', log = TRUE, sort = T, cols = c('gray','red','white','gray','darkseagreen', 'yellow', 'seagreen'))

p2 = DimPlot(data, group.by = 'seurat_clusters', reduction = 'umap.rna_0.5', label = T, cols = c('mediumaquamarine','seagreen','mediumseagreen','palegreen3','darkseagreen', 'cornflowerblue', 'royalblue', 'lightblue', 'palevioletred'))
VlnPlot(data[, data$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'seurat_clusters', log = TRUE, sort = T, cols = c('seagreen','mediumaquamarine','mediumseagreen','darkseagreen','palevioletred', 'royalblue', 'palegreen3', 'cornflowerblue', 'lightblue'))

p1 + p2 

VlnPlot(data, features = 'nFeature_RNA', group.by = 'seurat_clusters', log = TRUE, sort = T, cols = c('palevioletred', 'lightblue', 'royalblue', 'palegreen3', 'cornflowerblue', 'darkseagreen', 'mediumseagreen', 'mediumaquamarine', 'seagreen'))

FeaturePlot(data, features = 'ZAP70', order  = TRUE, cols = c("grey", "red"), pt.size = 0.5, reduction = 'umap.rna_0.5')


## barplot to observe the proportion of RT_annotations per cluster
proportions = as.data.frame((table(data$seurat_clusters, data$matchScore_RT_annotation) / rowSums(table(data$seurat_clusters, data$matchScore_RT_annotation)))*100) 

ggplot(proportions, aes(x = as.factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +  # Change position to "fill"
  labs(x = "Clusters", y = "Proportion", fill = "matchScore_RT_annotation") +
  theme_bw() +
  theme(legend.position = "bottom") +   
  scale_fill_manual(values = c("royalblue","cornflowerblue", "lightblue","darkseagreen","palegreen4", "palegoldenrod", "orange","palevioletred",'pink', "gray"))

table(data$matchScore_RT_annotation)

# CCND2hi RT       CCND2lo RT          CLL    CXCR4hiCD27lo    CXCR4loCD27hi 
#    10                7               17            11028             3865 

# CXCR4loCD27hi RT    MIR155HGhi RT RT proliferative    TP53INP1hi RT     unclassified 
#         1             1640             3624               97              418 


#ZAP70 proportions per sample

library(ggplot2)

proportions <- as.data.frame((table(data$seurat_clusters, data$ZAP_expression) / rowSums(table(data$seurat_clusters, data$ZAP_expression)))*100)
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

####
ZAP70_expression_plot

FeaturePlot(data, features = c("ZAP70", "CXCR4"), blend = TRUE, order = T)
FeaturePlot(data, features = c("ZAP70", "CD27"), blend = TRUE, order = T)
FeaturePlot(data, features = c("ZAP70", "TP53INP1"), blend = TRUE, order = T)

FeaturePlot(data, features = c("ZAP70", "CD38"), blend = TRUE, order = T)
# pacientes que expresan tp53 tienen mayor probabilidad de presentar una enfermedad que progresa más 
# rápidamente, requiere terapia, no responde bien a los tratamientos tradicionales

FeaturePlot(data, features = c("ZAP70", "TP53"), blend = TRUE, order = T)

########
# INTEGRATIONS (https://satijalab.org/seurat/articles/seurat5_integration)
data = SetIdent(data, value = data$orig.ident)

# 1. Anchor-based CCA integration (method=CCAIntegration)
data <- IntegrateLayers(object = data, method = CCAIntegration,
                        orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

data <- FindNeighbors(data, reduction = "integrated.cca", dims = 1:30)
data <- FindClusters(data, resolution = 0.5, cluster.name = "cca_clusters")
data <- RunUMAP(data, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
p1 = DimPlot(data,reduction = "umap.cca",group.by = "cca_clusters", combine = FALSE)

cca = DimPlot(data, reduction = "umap.cca", group.by = "RNA_snn_res.0.5", combine = FALSE, cols = c('mediumaquamarine','seagreen','mediumseagreen','palegreen3','darkseagreen', 'cornflowerblue', 'royalblue', 'lightblue', 'palevioletred'))

# Samples
DimPlot(data, reduction = "umap.cca", group.by = "orig.ident", combine = FALSE, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))
# RT annotation
DimPlot(data, reduction = "umap.cca", group.by = "matchScore_RT_annotation", combine = FALSE, cols = c('gray', 'gray', 'gray', 'darkgreen', 'palegreen3', 'gray', 'red', 'yellow', 'black', 'gray'))
# ZAP70 expression
FeaturePlot(data, features = 'ZAP70', reduction = 'umap.cca')


# 2. Anchor-based RPCA integration (method=RPCAIntegration)
data <- IntegrateLayers(object = data, method = RPCAIntegration,
                        orig.reduction = "pca", new.reduction = "integrated.rpca",
                        verbose = FALSE)

data <- FindNeighbors(data, reduction = "integrated.rpca", dims = 1:30)
data <- FindClusters(data, resolution = 0.5, cluster.name = "rpca_clusters")
data <- RunUMAP(data, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
p2 = DimPlot(data,reduction = "umap.rpca",group.by = "rpca_clusters",combine = FALSE)

rpca = DimPlot(data, reduction = "umap.rpca", group.by = "RNA_snn_res.0.5", combine = FALSE, cols = c('mediumaquamarine','seagreen','mediumseagreen','palegreen3','darkseagreen', 'cornflowerblue', 'royalblue', 'lightblue', 'palevioletred'))

# Samples
DimPlot(data, reduction = "umap.rpca", group.by = "orig.ident", combine = FALSE, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))
# RT annotation
DimPlot(data, reduction = "umap.rpca", group.by = "matchScore_RT_annotation", combine = FALSE, cols = c('gray', 'gray', 'gray', 'darkgreen', 'palegreen3', 'gray', 'red', 'yellow', 'black', 'gray'))
# ZAP70 expression
FeaturePlot(data, features = 'ZAP70', reduction = 'umap.rpca')

# 3. Harmony (method=HarmonyIntegration)
data <- IntegrateLayers(object = data, method = HarmonyIntegration,
                        orig.reduction = "pca", new.reduction = "harmonyoriginal",
                        verbose = FALSE)

data <- FindNeighbors(data, reduction = "harmonyoriginal", dims = 1:30)
data <- FindClusters(data, resolution = 0.5, cluster.name = "harmony_original_clusters")
data <- RunUMAP(data, reduction = "harmonyoriginal", dims = 1:30, reduction.name = "umap.harmony.original")
p3 = DimPlot(data,reduction = "umap.harmony.original",group.by = "harmony_original_clusters",combine = FALSE)

harmony = DimPlot(data, reduction = "umap.harmony.original", group.by = "RNA_snn_res.0.5", combine = FALSE, cols = c('mediumaquamarine','seagreen','mediumseagreen','palegreen3','darkseagreen', 'cornflowerblue', 'royalblue', 'lightblue', 'palevioletred'))

# Samples
DimPlot(data, reduction = "umap.harmony.original", group.by = "orig.ident", combine = FALSE, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))
# RT annotation
DimPlot(data, reduction = "umap.harmony.original", group.by = "matchScore_RT_annotation", combine = FALSE, cols = c('gray', 'gray', 'gray', 'darkgreen', 'palegreen3', 'gray', 'red', 'yellow', 'black', 'gray'))
# ZAP70 expression
FeaturePlot(data, features = 'ZAP70', reduction = 'umap.harmony.original')

# 4.FastMNN (method= FastMNNIntegration)
data <- IntegrateLayers(object = data, method = FastMNNIntegration,
                        new.reduction = "integrated.mnn",
                        verbose = FALSE)

data <- FindNeighbors(data, reduction = "integrated.mnn", dims = 1:30)
data <- FindClusters(data, resolution = 0.5, cluster.name = "mnn_clusters")
data <- RunUMAP(data, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")
p4 = DimPlot(data,reduction = "umap.mnn",group.by = "mnn_clusters",combine = FALSE)

mnn = DimPlot(data, reduction = "umap.mnn", group.by = "RNA_snn_res.0.5", combine = FALSE, cols = c('mediumaquamarine','seagreen','mediumseagreen','palegreen3','darkseagreen', 'cornflowerblue', 'royalblue', 'lightblue', 'palevioletred' p4))

# Samples
DimPlot(data, reduction = "umap.mnn", group.by = "orig.ident", combine = FALSE, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))
# RT annotation
DimPlot(data, reduction = "umap.mnn", group.by = "matchScore_RT_annotation", combine = FALSE, cols = c('gray', 'gray', 'gray', 'darkgreen', 'palegreen3', 'gray', 'red', 'yellow', 'black', 'gray'))
# ZAP70 expression
FeaturePlot(data, features = 'ZAP70', reduction = 'umap.mnn')


## All integration methods with original clusters

cca = DimPlot(data, reduction = "umap.cca", group.by = "RNA_snn_res.0.5", combine = FALSE, cols = c('mediumaquamarine','seagreen','mediumseagreen','palegreen3','darkseagreen', 'cornflowerblue', 'royalblue', 'lightblue', 'palevioletred'))
cca
rpca = DimPlot(data, reduction = "umap.rpca", group.by = "RNA_snn_res.0.5", combine = FALSE, cols = c('mediumaquamarine','seagreen','mediumseagreen','palegreen3','darkseagreen', 'cornflowerblue', 'royalblue', 'lightblue', 'palevioletred'))
rpca
harmony = DimPlot(data, reduction = "umap.harmony", group.by = "RNA_snn_res.0.5", combine = FALSE, cols = c('mediumaquamarine','seagreen','mediumseagreen','palegreen3','darkseagreen', 'cornflowerblue', 'royalblue', 'lightblue', 'palevioletred'))
harmony
mnn = DimPlot(data, reduction = "umap.mnn", group.by = "RNA_snn_res.0.5", combine = FALSE, cols = c('mediumaquamarine','seagreen','mediumseagreen','palegreen3','darkseagreen', 'cornflowerblue', 'royalblue', 'lightblue', 'palevioletred'))
mnn


## run harmony normal
data <- RunHarmony(data, "orig.ident")
data <- RunUMAP(data, dims = 1:20, reduction = 'harmony', reduction.name = 'umap.harmony')

DimPlot(data, reduction = "umap.harmony", group.by = "RNA_snn_res.0.5", combine = FALSE, cols = c('mediumaquamarine','seagreen','mediumseagreen','palegreen3','darkseagreen', 'cornflowerblue', 'royalblue', 'lightblue', 'palevioletred'))

data <- FindNeighbors(data, reduction = "harmony", dims = 1:30)
data <- FindClusters(data, resolution = 0.5, cluster.name = "harmony_clusters")
data <- RunUMAP(data, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

DimPlot(data,reduction = "umap.harmony",group.by = "harmony_clusters",combine = FALSE)
DimPlot(data, reduction = "umap.harmony", group.by = "matchScore_RT_annotation", combine = FALSE, cols = c('gray','gray','gray','seagreen','palegreen3', 'gray', 'red', 'yellow', 'black', 'gray'))
DimPlot(data,reduction = "umap.harmony",group.by = "orig.ident",combine = FALSE)
DimPlot(data, reduction = "umap.harmony", group.by = "RNA_snn_res.0.5", combine = FALSE, cols = c('mediumaquamarine','seagreen','mediumseagreen','palegreen3','darkseagreen', 'cornflowerblue', 'royalblue', 'lightblue', 'palevioletred', 'pink'))

FeaturePlot(data, features = 'ZAP70', reduction = 'umap.harmony')

## barplot to observe the proportion of RT_annotations per cluster
proportions = as.data.frame((table(data$harmony_clusters, data$matchScore_RT_annotation) / rowSums(table(data$harmony_clusters, data$matchScore_RT_annotation)))*100) 

ggplot(proportions, aes(x = as.factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +  # Change position to "fill"
  labs(x = "Clusters", y = "Proportion", fill = "matchScore_RT_annotation") +
  theme_bw() +
  theme(legend.position = "bottom") +   
  scale_fill_manual(values = c("royalblue","cornflowerblue", "lightblue","darkseagreen","palegreen4", "palegoldenrod", "orange","palevioletred",'pink', "gray"))


'''
# 5. scVI (method=scVIIntegration) !!!!!!!
library(reticulate)
data <- IntegrateLayers(
  object = data, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/opt/homebrew/Caskroom/mambaforge/base/envs/scvi-env/", verbose = FALSE
)
'''

## 7.A TRYING OTHER RESOLUTIONS (1)

data <- FindClusters(data, resolution = 1) #Partition the KNN graph into clusters.
table(data@active.ident)
#   0    1    2    3    4    5    6    7    8    9   10   11   12 
# 4541 4051 3886 1960 1722 1418  911  735  550  338  311  165  119

data <- RunUMAP(data, dims = 1:20, reduction.name = 'umap.rna.1')

DimPlot(data, group.by = 'seurat_clusters', reduction = 'umap.harmony', label = T)
DimPlot(data, group.by = 'matchScore_RT_annotation', reduction = 'umap.harmony', label = F, cols = c('gray','gray','gray','darkseagreen','seagreen', 'gray', 'red', 'yellow', 'black','gray'))

p1 | p2

table(data$RNA_snn_res.0.8)

## 7.B TRYING OTHER RESOLUTIONS (1.2)

data <- FindClusters(data, resolution = 1.2) #Partition the KNN graph into clusters.
table(data@active.ident)
#   0    1    2    3    4    5    6    7    8    9   10   11   12    13   14   
# 4443 3862 2366 1901 1788 1637 1411  920  763  552  338  314  162  131  119 

data <- RunUMAP(data, dims = 1:20, reduction.name = 'umap.rna.1.2')

p1 = DimPlot(data, group.by = 'RNA_snn_res.1.2', reduction = 'umap.harmony', label = T)
p2 = DimPlot(data, group.by = 'matchScore_RT_annotation', reduction = 'umap.harmony', label = F, cols = c('gray','gray','gray','darkseagreen','seagreen', 'gray', 'red', 'yellow', 'black','gray'))
p3 = DimPlot(data, group.by = 'orig.ident', reduction = 'umap.harmony', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen'))

p1 | p2 | p3

'''
model = readRDS('after_meeting/model_RT_annotations.rds')

data <- JoinLayers(data)
out <- identity_map(scale.data = data@assays$RNA$scale.data ,model = model,top)
ids <- out$ids 

table(ids)

# CCND2lo RT      CLL    CXCR4hiCD27lo    CXCR4loCD27hi 
#       2          5         2943             1368 

# MIR155HGhi RT    RT proliferative    TP53INP1hi RT     unclassified 
#      664              852               25                 108 

data$matchScore_RT_annotation_1 <- ids # add to my metadata the annotation of my model trained with the RT dataset
'''

# unify clusters nearby
library(plyr)
data$my_clusters <- revalue(data$RNA_snn_res.1.2, c("0"="0", "1"="0", "2"= "1", "3"= "1", "4" = "1","5" = "2","6" = "3","7" = "3","8" = "2","9" = "4", "10" = "5","11" = "6", "12" = "7", "13" = "0","14" = "8"))
                          
# reannotation 
p1 = DimPlot(data, group.by = 'my_clusters', reduction = 'umap.harmony', label = T)
# matchscore
p2 = DimPlot(data, group.by = 'matchScore_RT_annotation', reduction = 'umap.harmony', label = F, cols = c('gray','gray','gray','darkseagreen','seagreen', 'gray', 'red', 'yellow', 'black','gray'))
# samples
p3 = DimPlot(data, group.by = 'orig.ident', reduction = 'umap.harmony', label = T)


## barplot to observe the proportion of RT_annotations per cluster
proportions = as.data.frame((table(data$my_clusters, data$matchScore_RT_annotation_1) / rowSums(table(data$my_clusters, data$matchScore_RT_annotation_1)))*100) 

p4 = ggplot(proportions, aes(x = as.factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +  # Change position to "fill"
  labs(x = "Clusters", y = "Proportion", fill = "matchScore_RT_annotation") +
  theme_bw() +
  theme(legend.position = "bottom") +   
  scale_fill_manual(values = c("gray","gray", "gray","darkseagreen","seagreen", "gray", "red","yellow",'black', "gray"))


p5 = DotPlot(data, features = c('CXCR4', 'CD24','CD27', 'MIR155HG','CCND2', 'PCNA', 'MKI67','MZB1','IGHM', 'XBP1'), group.by = 'my_clusters') + RotatedAxis()

p1 + p5

### Generate a new model with the RT dataset excluding the RT clusters

# Remove the richter clusters from my data
richter_clusters <- c("CCND2hi RT", "CCND2lo RT", "CXCR4loCD27hi RT", "MIR155HGhi RT", "RT", "RT proliferative", "RT quiescent", "TP53INP1hi RT") # identify CD4+ T cells and CD8+ Cytotoxic T cells
richter_clusters <- which(data_rt@meta.data$annotation_fina %in% richter_clusters) # identify rt (positions)
data_rt_no_clusters <- subset(data_rt, cells = -richter_clusters) # subset to remove T cells

Idents(data_rt_no_clusters) = data_rt_no_clusters$annotation_fina
markers <- FindAllMarkers(data_rt_no_clusters, only.pos = T)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=100)

scaled_ref <- ScaleData(data_rt_no_clusters, features = unlist(top))

## Generate the model 
model <- train_model(scale.data = scaled_ref,clus = data_rt_no_clusters$annotation_final, gene_cl.ref = top,prop = 0.75)

## Cell projection in my data
out <- identity_map(scale.data = data@assays$RNA$scale.data,model = model,top)
 
### cell identities
ids <- out$ids 

data$annotation_noRT_clusters = ids

DimPlot(data, group.by = 'annotation_noRT_clusters', reduction = 'umap.harmony', label = F, cols = c('grey', 'black', 'palegoldenrod', 'red', 'green4'))

table(data$annotation_noRT_clusters)

# CD83loMIR155HGhi     CLL    CXCR4hiCD27lo    CXCR4loCD27hi     unclassified 
#     2               93            18066             2520               26 

### Find all markers of my new clusters to annotate them (+ bibliography)


## 9. FINDING DIFFERENTIALLY EXPRESSED FEATURES (cluster biomarkers)
Idents(data) = data$my_clusters

markers <- FindAllMarkers(data, only.pos = T)
write.csv(markers, file = "after_meeting/markers_myclusters.csv")

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

DoHeatmap(data, features = top10$gene) + NoLegend()

top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=100)
write.csv(top, file = "after_meeting/top100markers_myclusters.csv")

top_table = as.data.frame(top)

top_table$X0 
# I consider important: IGKC, CCND3, MAST4, NOTCH2, GNG7.  (cluster 'MIR155HGhi' MIR155HG)

top_table$X1
# I consider important: IGLC2, IGLC3, IGHM, CXCR4, JUND, JUNB, CD69, XBP1 (cluster 'CD83hiMIR155HGhi' = IGLC3)

top_table$X2
# I consider important: CD69, IGLC3, IGLC2, CD83, JUND, JUN, MAP2?, XBP1

top_table$X3
# I consider important: MAST4, GNG7, NOTCH2

top_table$X4
# I consider important: IGKC, CCDC88A, MS4A1, CD74?

top_table$X5
# I consider important: DNAJB9, XBP1

top_table$X6
# I consider important: CD3E, GZMM, STAT4, GNLY, CD8A, TNFAIP3, ZAP70, IL7R, CD3D, NKG7

top_table$X7
# I consider important: STAT4, IL7R, CD3D, TNFAIP3, CD3E, GZMM, CD3G, 

top_table$X8
# I consider important: MIR3681HG


# The unmutated group contained significantly higher percentages of leukemic cells expressing CD38, CD69, and CD40
# LOOK FOR TCELL MARKERS
DotPlot(data, features = c('CXCR4', 'CD24','CD27', 'MIR155HG','CCND2', 'PCNA', 'MKI67','MZB1','IGHM', 'XBP1', 'CD8A'), group.by = 'my_clusters') + RotatedAxis()
FeaturePlot(data, features = c('CD3E', 'CD3G', 'CD8A', 'IL7R'), reduction = 'umap.harmony')

### GO TO noCD3analysis.R script to follow analysis
# Per sample (RT dataset)
patient19 = readRDS('RT/paciente19.rds')
patient12 = readRDS('RT/paciente12.rds')
patient365 = readRDS('RT/paciente365.rds')
patient3299 = readRDS('RT/paciente3299.rds')



  