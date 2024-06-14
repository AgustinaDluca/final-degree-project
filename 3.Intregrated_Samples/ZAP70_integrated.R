---
  title: "all samples"
author: "Agustina De Luca"
date: "2024-02-20"
output: html_document
---
  
### ZAP70 ###
  
FeaturePlot(data, features = 'ZAP70',order = T) # location in the clusters
VlnPlot(data, features = 'ZAP70', slot = "counts", log = TRUE, cols = c('#118AB2','#06D6A0','#8338EC','orange','#FF6B6B', 'violet', '#F4D35E', 'lightblue', 'darkgreen', 'deepskyblue2')) # RNA counts

## Label cells that express the ZAP70 marker in more than one count per cell
# Create a column in the metadata for ZAP expression, here we will evaluate if a cell is ZAP positive or not
data$ZAP_expression <- ifelse(data$RNA$counts['ZAP70', ] >= 1, "ZAP+", "ZAP-")

table(data$ZAP_expression)
# ZAP-  ZAP+ 
# 21230  1146 

table(data$seurat_clusters, data$ZAP_expression)
#  ZAP- ZAP+
#0 5728   35
#1 5458  275
#2 4927    1
#3 1077  247
#4 1034  214
#5  875  213
#6 1005   15
#7  584    4
#8  309  131
#9  233   11

## PLOT PER CLUSTER ZAP+/ZAP-
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

ZAP70_expression_plot

data$ZAP_expression <- as.factor(data$ZAP_expression)

# EXPRESSION
FeaturePlot(data, features = 'ZAP70', order  = TRUE, cols = c("grey", "red"), pt.size = 0.5)

# COLOR CELLS THAT ARE ZAP+ 
DimPlot(data, group.by = "ZAP_expression", cols = c("lightgray", "red"), pt.size = 0.1, order = T)


##Create a violin plot for ZAP70 in ZAP70+ cells
#S_[, S_$ZAP_expression == 'ZAP+'] -> FILTER cells that are ZAP70+                                                                                                                                    
VlnPlot(data[, data$ZAP_expression == 'ZAP+', drop = FALSE], features = 'ZAP70', group.by = 'cell_type', log = TRUE, sort = T, cols = c('pink','palevioletred','darkgreen','mediumaquamarine','mediumseagreen','#F4D35E', 'darkseagreen', 'salmon1'))

#ZAP+/- clusters
data$cluster_expression <- paste(data$seurat_clusters, data$ZAP_expression, sep = "_")

##Dot plot of markers per cluster of each sample (agarre los primeros dos markers de cada clsuter)
DotPlot(data, features = c(...), col.min = 0.8, dot.scale = 5, scale.by = 'size') + RotatedAxis()


##Dot plot of clusters/ZAP70 with common markers
Idents(data) = as.factor(data$sub.cluster)
DotPlot(data, features = c("IL7R", "CCR7", "S100A4", "CD4","CD8A", "CD3E", "CD3G", "GNLY", "NKG7", "MS4A1", "FCGR3A", "MS4A7", "CD14", "LYZ", "FCER1A", "CST3", "PPBP"), col.min = 0) + RotatedAxis()

## Heatmap w/ DE-genes from ZAP+ and ZAP-
## Get top 100 from ZAP+/ZAP-

data <- SetIdent(data, value = data$ZAP_expression)
DE_genes <- FindMarkers(object = data, ident.1 = 'ZAP+', ident.2 = 'ZAP-') # Run differential expression analysis
#write.csv(DE_genes, file = "DE_genes_integrated.csv", row.names = T)

integration_ZAP70_pos <- head(DE_genes$X, 500)
integration_ZAP70_neg <- tail(DE_genes$X, 500)

integration_ZAP70 <- data.frame(integration_ZAP70_pos, integration_ZAP70_neg)
names(integration_ZAP70) <- c("ZAP+", "ZAP-")

DE_genes = read.csv('DE_genes_integrated.csv')

head(integration_ZAP70, 10)

## Coexpression of ZAP70+ and ZAP70- with the common markers

#B-cells
FeaturePlot(data, features = c("ZAP70", "MS4A1"), blend = TRUE, order = T)
#NK cells
FeaturePlot(data, features = c("ZAP70", "GNLY"), blend = TRUE, order = T)

#T-celL
FeaturePlot(data, features = c("ZAP70", "CD3D"), blend = TRUE, order = T)#
FeaturePlot(data, features = c("ZAP70", "CD3E"), blend = TRUE, order = T)#
FeaturePlot(data, features = c("ZAP70", "CD3G"), blend = TRUE, order = T)#
FeaturePlot(data, features = c("ZAP70", "CD8A"), blend = TRUE, order = T)
FeaturePlot(data, features = c("ZAP70", "CD4"), blend = TRUE, order = T)
FeaturePlot(data, features = c("ZAP70", "IL7R"), blend = TRUE, order = T)

#B-cells
FeaturePlot(data, features = c("ZAP70", "CD79A"), blend = TRUE, order = T)
FeaturePlot(data, features = c("ZAP70", "CD79B"), blend = TRUE, order = T)

#Monocytes
FeaturePlot(data, features = c("ZAP70", "FCGR3A"), blend = TRUE, order = T)#

#Langerhans cells
FeaturePlot(data, features = c("ZAP70", "LYZ"), blend = TRUE, order = T)

# Notably, ZAP-70 level in ALL is associated with CD38 expression, but no correlation was observed to specific cytogenetic abnormalities 
FeaturePlot(data, features = c("ZAP70", "CD38"), blend = TRUE, order = T)

FeaturePlot(data, features = c("ZAP70", "STAT4"), blend = TRUE, order = T)
FeaturePlot(data, features = c("ZAP70", "THEMIS"), blend = TRUE, order = T)

## DE-genes from ZAP+ and ZAP- per cluster!!! IMPORTANT!!!
## Get top 100 from ZAP+/ZAP-

data <- SetIdent(data, value = as.factor(data$cluster_expression))

DE_genes_cluster0 <- FindMarkers(object = data, ident.1 = '0_ZAP+', ident.2 = '0_ZAP-') # Run differential expression analysis
DE_genes_cluster0 = DE_genes_cluster0[order(-DE_genes_cluster0[,2], DE_genes_cluster0[,1] ),]
write.csv(DE_genes_cluster0, file = "Documents/Documents/INTEGRATED DATA/DE_genes_cluster0.csv", row.names = T)

DE_genes_cluster0 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_cluster0_pos # aca tengo todos los ZAP70+ genes
write.csv(DE_cluster0_pos, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_cluster0_pos.csv", row.names = T)

DE_genes_cluster0 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC < -0.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_cluster0_neg # aca tengo todos los ZAP70- genes
write.csv(DE_cluster0_neg, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_cluster0_neg.csv", row.names = T)


DE_genes_cluster3 <- FindMarkers(object = data, ident.1 = '3_ZAP+', ident.2 = '3_ZAP-') # Run differential expression analysis
DE_genes_cluster3 = DE_genes_cluster3[order(-DE_genes_cluster3[,2], DE_genes_cluster3[,1] ),]
write.csv(DE_genes_cluster3, file = "Documents/Documents/INTEGRATED DATA/DE_genes_cluster3.csv", row.names = T)

DE_genes_cluster3 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_cluster3_pos # aca tengo todos los ZAP70+ genes
write.csv(DE_cluster3_pos, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_cluster3_pos.csv", row.names = T)

DE_genes_cluster3 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC < -0.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_cluster3_neg # aca tengo todos los ZAP70- genes
write.csv(DE_cluster3_neg, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_cluster3_neg.csv", row.names = T)



DE_genes_cluster4 <- FindMarkers(object = data, ident.1 = '4_ZAP+', ident.2 = '4_ZAP-') # Run differential expression analysis
DE_genes_cluster4 = DE_genes_cluster4[order(-DE_genes_cluster4[,2], DE_genes_cluster4[,1] ),]
write.csv(DE_genes_cluster4, file = "Documents/Documents/INTEGRATED DATA/DE_genes_cluster4.csv", row.names = T)

DE_genes_cluster4 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_cluster4_pos # aca tengo todos los ZAP70+ genes
write.csv(DE_cluster4_pos, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_cluster4_pos.csv", row.names = T)

DE_genes_cluster4 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC < -0.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_cluster4_neg # aca tengo todos los ZAP70- genes
write.csv(DE_cluster4_neg, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_cluster4_neg.csv", row.names = T)



DE_genes_cluster6 <- FindMarkers(object = data, ident.1 = '6_ZAP+', ident.2 = '6_ZAP-') # Run differential expression analysis
DE_genes_cluster6 = DE_genes_cluster6[order(-DE_genes_cluster6[,2], DE_genes_cluster6[,1] ),]
write.csv(DE_genes_cluster6, file = "Documents/Documents/INTEGRATED DATA/DE_genes_cluster6.csv", row.names = T)

DE_genes_cluster6 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_cluster6_pos # aca tengo todos los ZAP70+ genes
write.csv(DE_cluster6_pos, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_cluster6_pos.csv", row.names = T)

DE_genes_cluster6 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC < -0.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_cluster6_neg # aca tengo todos los ZAP70- genes
write.csv(DE_cluster6_neg, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_cluster6_neg.csv", row.names = T)



DE_genes_cluster8 <- FindMarkers(object = data, ident.1 = '8_ZAP+', ident.2 = '8_ZAP-') # Run differential expression analysis
DE_genes_cluster8 = DE_genes_cluster8[order(-DE_genes_cluster8[,2], DE_genes_cluster8[,1] ),]
write.csv(DE_genes_cluster8, file = "Documents/Documents/INTEGRATED DATA/DE_genes_cluster8.csv", row.names = T)
DE_genes_cluster8 = read.csv("Documents/Documents/INTEGRATED DATA/DE_genes_cluster8.csv")

DE_genes_cluster8 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_cluster8_pos # aca tengo todos los ZAP70+ genes
write.csv(DE_cluster8_pos, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_cluster8_pos.csv", row.names = T)

DE_genes_cluster8 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC < -0.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_cluster8_neg # aca tengo todos los ZAP70- genes
write.csv(DE_cluster8_neg, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_cluster8_neg.csv", row.names = T)



DE_genes_cluster9 <- FindMarkers(object = data, ident.1 = '9_ZAP+', ident.2 = '9_ZAP-') # Run differential expression analysis
DE_genes_cluster9 = DE_genes_cluster9[order(-DE_genes_cluster9[,2], DE_genes_cluster9[,1] ),]
write.csv(DE_genes_cluster9, file = "Documents/Documents/INTEGRATED DATA/DE_genes_cluster9.csv", row.names = T)
DE_genes_cluster9 = read.csv("Documents/Documents/INTEGRATED DATA/DE_genes_cluster9.csv")

DE_genes_cluster9 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_cluster9_pos # aca tengo todos los ZAP70+ genes
write.csv(DE_cluster9_pos, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_cluster9_pos.csv", row.names = T)

DE_genes_cluster9 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC < -0.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_cluster9_neg # aca tengo todos los ZAP70- genes
write.csv(DE_cluster9_neg, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_cluster9_neg.csv", row.names = T)

DE_cluster0_pos <- read.csv('Documents/Documents/scRNA   intregrated data/DE genes/DE_cluster0_pos.csv')
DE_cluster0_neg <- read.csv('Documents/Documents/scRNA   intregrated data/DE genes/DE_cluster0_neg.csv')

DE_cluster3_pos <- read.csv('Documents/Documents/scRNA   intregrated data/DE genes/DE_cluster3_pos.csv')
DE_cluster3_neg <- read.csv('Documents/Documents/scRNA   intregrated data/DE genes/DE_cluster3_neg.csv')

DE_cluster4_pos <- read.csv('Documents/Documents/scRNA   intregrated data/DE genes/DE_cluster4_pos.csv')
DE_cluster4_neg <- read.csv('Documents/Documents/scRNA   intregrated data/DE genes/DE_cluster4_neg.csv')

DE_cluster6_pos <- read.csv('Documents/Documents/scRNA   intregrated data/DE genes/DE_cluster6_pos.csv')
DE_cluster6_neg <- read.csv('Documents/Documents/scRNA   intregrated data/DE genes/DE_cluster6_neg.csv')

DE_cluster8_pos <- read.csv('Documents/Documents/scRNA   intregrated data/DE genes/DE_cluster8_pos.csv')
DE_cluster8_neg <- read.csv('Documents/Documents/scRNA   intregrated data/DE genes/DE_cluster8_neg.csv')

DE_cluster9_pos <- read.csv('Documents/Documents/scRNA   intregrated data/DE genes/DE_cluster9_pos.csv')
DE_cluster9_neg <- read.csv('Documents/Documents/scRNA   intregrated data/DE genes/DE_cluster9_neg.csv')

### GO ENTRICHMENT ANALYSIS
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html

# -- CLUSTER0 -- #

# ZAP70+ 

# Entrez IDs
ZAP_entrez_pos <- bitr(DE_cluster0_pos$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# muchos no mapean justamente porque son unknown markers por eso quedan tan pocos

# GO analysis
enrichment_results_cluster0_pos <- enrichGO(gene = ZAP_entrez_pos$ENTREZID,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                            pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.05)
# Visualization
barplot(enrichment_results_cluster0_pos)
dotplot(enrichment_results_cluster0_pos)


# ZAP70-

# Entrez IDs
ZAP_entrez_neg <- bitr(DE_cluster0_neg$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_cluster0_neg <- enrichGO(gene = ZAP_entrez_neg$ENTREZID,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                            pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                            pvalueCutoff = 0.08,
                                            qvalueCutoff = 0.08)
# Visualization
barplot(enrichment_results_cluster0_neg)
dotplot(enrichment_results_cluster0_neg)



# -- CLUSTER3 -- #

# ZAP70+ 

# Entrez IDs
ZAP_entrez_pos <- bitr(DE_cluster3_pos$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_cluster3_pos <- enrichGO(gene = ZAP_entrez_pos$ENTREZID,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                            pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                            pvalueCutoff = 0.08,
                                            qvalueCutoff = 0.08)
# Visualization
barplot(enrichment_results_cluster3_pos)
dotplot(enrichment_results_cluster3_pos)


# ZAP70-

# Entrez IDs
ZAP_entrez_neg <- bitr(DE_cluster3_neg$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_cluster3_neg <- enrichGO(gene = ZAP_entrez_neg$ENTREZID,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                            pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.05)
# Visualization
barplot(enrichment_results_cluster3_neg)
dotplot(enrichment_results_cluster3_neg)

# -- CLUSTER4 -- #

# ZAP70+ 

# Entrez IDs
ZAP_entrez_pos <- bitr(DE_cluster4_pos$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_cluster4_pos <- enrichGO(gene = ZAP_entrez_pos$ENTREZID,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                            pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                            pvalueCutoff = 0.08,
                                            qvalueCutoff = 0.08)
# Visualization
barplot(enrichment_results_cluster4_pos)
dotplot(enrichment_results_cluster4_pos)


# ZAP70-

# Entrez IDs
ZAP_entrez_neg <- bitr(DE_cluster4_neg$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_cluster4_neg <- enrichGO(gene = ZAP_entrez_neg$ENTREZID,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                            pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                            pvalueCutoff = 0.3,
                                            qvalueCutoff = 0.3)
# Visualization
barplot(enrichment_results_cluster4_neg)
dotplot(enrichment_results_cluster4_neg)

# -- CLUSTER6 -- #

# ZAP70+ 

# Entrez IDs
ZAP_entrez_pos <- bitr(DE_cluster6_pos$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_cluster6_pos <- enrichGO(gene = ZAP_entrez_pos$ENTREZID,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                            pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.05)
# Visualization
barplot(enrichment_results_cluster6_pos)
dotplot(enrichment_results_cluster6_pos)


# ZAP70-

# Entrez IDs
ZAP_entrez_neg <- bitr(DE_cluster6_neg$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_cluster6_neg <- enrichGO(gene = ZAP_entrez_neg$ENTREZID,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                            pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                            pvalueCutoff = 0.1,
                                            qvalueCutoff = 0.1)
# Visualization
barplot(enrichment_results_cluster6_neg)
dotplot(enrichment_results_cluster6_neg)

# -- CLUSTER8 -- #

# ZAP70+ 

# Entrez IDs
ZAP_entrez_pos <- bitr(DE_cluster8_pos$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_cluster8_pos <- enrichGO(gene = ZAP_entrez_pos$ENTREZID,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                            pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                            pvalueCutoff = 0.08,
                                            qvalueCutoff = 0.08)
# Visualization
barplot(enrichment_results_cluster8_pos)
dotplot(enrichment_results_cluster8_pos)


# ZAP70-

# Entrez IDs
ZAP_entrez_neg <- bitr(DE_cluster8_neg$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_cluster8_neg <- enrichGO(gene = ZAP_entrez_neg$ENTREZID,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                            pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.05)
# Visualization
barplot(enrichment_results_cluster8_neg)
dotplot(enrichment_results_cluster8_neg)

# -- CLUSTER9 -- #

# ZAP70+ 

# Entrez IDs
ZAP_entrez_pos <- bitr(DE_cluster9_pos$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_cluster9_pos <- enrichGO(gene = ZAP_entrez_pos$ENTREZID,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                            pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                            pvalueCutoff = 0.3,
                                            qvalueCutoff = 0.3)
# Visualization
barplot(enrichment_results_cluster9_pos)
dotplot(enrichment_results_cluster9_pos)


# ZAP70-

# Entrez IDs
ZAP_entrez_neg <- bitr(DE_cluster9_neg$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_cluster9_neg <- enrichGO(gene = ZAP_entrez_neg$ENTREZID,
                                            OrgDb = org.Hs.eg.db,
                                            ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                            pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                            pvalueCutoff = 0.3,
                                            qvalueCutoff = 0.3)
# Visualization
barplot(enrichment_results_cluster9_neg)
dotplot(enrichment_results_cluster9_neg)

zero = read.csv('Documents/Documents/INTEGRATED DATA/DE genes antes/DE_genes_cluster0.csv')

zero_pos = head(zero, 500)
zero_neg = tail(zero, 500)

zero_pos %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 1.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> zero_pos_filter # aca tengo todos los ZAP70+ genes

zero_neg %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC < -1.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> zero_neg_filter # aca tengo todos los ZAP70+ genes

three = read.csv('Documents/Documents/INTEGRATED DATA/DE genes antes/DE_genes_cluster3.csv')

three_pos = head(three, 500)
three_neg = tail(three, 500)

three_pos %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 1.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> three_pos_filter # aca tengo todos los ZAP70+ genes

three_neg %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC < -1.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> three_neg_filter # aca tengo todos los ZAP70+ genes

four = read.csv('Documents/Documents/INTEGRATED DATA/DE genes antes/DE_genes_cluster4.csv')


six = read.csv('Documents/Documents/INTEGRATED DATA/DE genes antes/DE_genes_cluster6.csv')
eight = read.csv('Documents/Documents/INTEGRATED DATA/DE genes antes/DE_genes_cluster8.csv')
nine = read.csv('Documents/Documents/INTEGRATED DATA/DE genes antes/DE_genes_cluster9.csv')


ZAP_entrez_pos <- bitr(zero_pos_filter$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO analysis
agusss <- enrichGO(gene = ZAP_entrez_pos$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                   pAdjustMethod = "BH", # Benjamini & Hochberg correction
                   pvalueCutoff = 0.08,
                   qvalueCutoff = 0.08)
# Visualization
barplot(enrichment_results_cluster0_pos)

# Volcano Plot
library(ggplot2)


variable <- agus %>% 
  mutate(significant1= ifelse(abs(avg_log2FC) > 1.5 & p_val < 0.05, "Significant", "Not significant"))

ggplot(variable, aes(x=avg_log2FC, y=-log(p_val), color = significant1)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("red", "grey50")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Log2 Fold Change", y = "P-value") +
  ggtitle("Volcano plot of differentially expressed genes") +
  theme(plot.title = element_text(hjust = 0.5))


data <- SetIdent(data, value = as.factor(data$ZAP_expression))

DE_genes_sample <- FindMarkers(object = data, ident.1 = 'ZAP+', ident.2 = 'ZAP-') # Run differential expression analysis
DE_genes_sample = DE_genes_sample[order(-DE_genes_sample[,2], DE_genes_sample[,1] ),]
write.csv(DE_genes_sample, file = "Documents/Documents/INTEGRATED DATA/DE_genes_sample.csv", row.names = T)

DE_genes_sample %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 1.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_genes_sample_pos # aca tengo todos los ZAP70+ genes
write.csv(DE_genes_sample_pos, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_genes_sample_pos.csv", row.names = T)

DE_genes_sample %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC < -1.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_genes_sample_neg # aca tengo todos los ZAP70- genes
write.csv(DE_genes_sample_neg, file = "Documents/Documents/INTEGRATED DATA/DE genes/DE_genes_sample_neg.csv", row.names = T)

DE_genes_sample_pos = read.csv('Documents/Documents/MULTIOME/scRNA   intregrated data/DE genes/DE_genes_sample_pos.csv')
DE_genes_sample_neg = read.csv('Documents/Documents/MULTIOME/scRNA   intregrated data/DE genes/DE_genes_sample_neg.csv')
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
                                          pvalueCutoff = 0.2,
                                          qvalueCutoff = 0.2)
# Visualization
barplot(enrichment_results_sample_neg)
dotplot(enrichment_results_sample_neg)


markers_rna_zap70 <- presto:::wilcoxauc.Seurat(X = data, group_by = 'ZAP_expression', assay = 'data', seurat_assay = 'RNA') # get RNA markers
markers_motifs_zap70 <- presto:::wilcoxauc.Seurat(X = data, group_by = 'ZAP_expression', assay = 'data', seurat_assay = 'chromvar') # get motifs

motif.names_zap70 <- markers_motifs_zap70$feature # extract the motifs names

colnames(markers_rna_zap70) <- paste0("RNA.", colnames(markers_rna_zap70)) # state the names of the columns in markers_RNA (variables)
colnames(markers_motifs_zap70) <- paste0("motif.", colnames(markers_motifs_zap70)) # state the names of the columns in markers_RNA (variables)

markers_rna_zap70$gene <- markers_rna_zap70$RNA.feature # the last column are the features, now genes from the RNA data
markers_motifs_zap70$gene <- ConvertMotifID(data, id = motif.names_zap70) # the last column are the motif names, now with the correspoinding ID

# a simple function to implement the procedure above for zap70 positive and negative
topTFs_zap70 <- function(celltype, padj.cutoff = 1e-2) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna_zap70, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
    arrange(-RNA.auc)
  ctmarkers_motif <- dplyr::filter(
    markers_motifs_zap70, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
    arrange(-motif.auc)
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(2, 11, 6, 7)], 
    y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
  )
  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, -avg_auc)
  return(top_tfs)
}

ZAP.neg = topTFs_zap70('ZAP-')

#  RNA.group  gene   RNA.auc     RNA.pval motif.group motif.feature motif.auc   motif.pval   avg_auc
# 1      ZAP- MEF2C 0.5850738 1.141005e-24        ZAP-      MA0497.1 0.6169454 1.064377e-40 0.6010096
# 2      ZAP-  TCF4 0.5428688 9.418675e-07        ZAP-      MA0830.2 0.5342061 9.339052e-05 0.5385374
# 3      ZAP-  MXI1 0.4730487 3.496565e-06        ZAP-      MA1108.2 0.5866100 4.467347e-23 0.5298294

ZAP.neg = topTFs_zap70('ZAP-')

#   RNA.group gene   RNA.auc      RNA.pval motif.group motif.feature motif.auc   motif.pval   avg_auc
# 1      ZAP+ JUND 0.6683968 4.045555e-105        ZAP+      MA0491.2 0.6550394 3.560626e-70 0.6617181
# 2      ZAP+  FOS 0.5957995 1.586899e-158        ZAP+      MA0476.1 0.6547398 6.536083e-70 0.6252696
# 3      ZAP+ JUNB 0.5955830  1.565369e-47        ZAP+      MA0490.2 0.6465879 6.280320e-63 0.6210854
# 4      ZAP+  REL 0.5633523  2.755121e-16        ZAP+      MA0101.1 0.6658743 4.691994e-80 0.6146133
# 5      ZAP+ RELA 0.5537619  1.590847e-50        ZAP+      MA0107.1 0.6653521 1.454595e-79 0.6095570


MotifPlot(object = data,
          motifs = ZAP.neg$motif.feature,
          assay = 'ATAC')
 

data@active.assay = 'RNA'
gene_plot <- FeaturePlot(data, features = c("MEF2C",'TCF4','MXI1') , reduction = 'wnn.umap')
data@active.assay = 'chromvar'
motif_plot <- FeaturePlot(data, features = c("MA0497.1",'MA0830.2', 'MA1108.2'), min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot


MotifPlot(object = data,
          motifs = ZAP.pos$motif.feature,
          assay = 'ATAC')


data@active.assay = 'RNA'
gene_plot <- FeaturePlot(data, features = c(ZAP.pos$gene[1:5]) , reduction = 'umap.rna')
data@active.assay = 'chromvar'
motif_plot <- FeaturePlot(data, features = c(ZAP.pos$motif.feature[1:5]), min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'umap.rna')
gene_plot | motif_plot


DefaultAssay(data) <- 'ATAC'

Idents(data) = as.factor(data$ZAP_expression)

da_peaks_ZAP70 <- FindMarkers(object = data, ident.1 = 'ZAP+', ident.2 = 'ZAP-') # Run differential expression analysis
head(da_peaks_ZAP70)

da_peaks_ZAP70 = da_peaks_ZAP70[order(-da_peaks_ZAP70[,2], da_peaks_ZAP70[,1] ),]

da_peaks_ZAP70 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 1.5) %>%
  slice_head(n = 100) %>% 
  ungroup() -> da_peaks_ZAP70_pos # aca tengo todos los ZAP70+ genes
write.csv(da_peaks_ZAP70_pos, file = "Documents/Documents/MULTIOME/da_peaks_ZAP70_pos.csv", row.names = T)

da_peaks_ZAP70 %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC < -1.5) %>%
  slice_head(n = 100) %>% 
  ungroup() -> da_peaks_ZAP70_neg # aca tengo todos los ZAP70- genes
write.csv(da_peaks_ZAP70_neg, file = "Documents/Documents/MULTIOME/da_peaks_ZAP70_neg.csv", row.names = T)

# Test enrichment
enriched.motifs_pos <- FindMotifs(object = data, features = rownames(da_peaks_ZAP70_pos[1:5]))

enriched.motifs_neg <- FindMotifs(object = data, features = rownames(da_peaks_ZAP70_neg[1:5]))

MotifPlot(object = data,motifs = head(rownames(enriched.motifs_pos)))
MotifPlot(object = data,motifs = head(rownames(enriched.motifs_neg)))

data@active.assay = 'ATAC'
motif.name <- ConvertMotifID(data, name = 'ZAP70')
data@active.assay = 'RNA'
gene_plot <- FeaturePlot(data, features = "TEAD3", reduction = 'wnn.umap')
data@active.assay = 'chromvar'
motif_plot <- FeaturePlot(data, features = 'MA0808.1', min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot
