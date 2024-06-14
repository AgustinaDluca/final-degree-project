#Richter tranformation patient

# Github: https://github.com/massonix/richter_transformation
# Zenodo: https://zenodo.org/records/6631966


# Paper interesante donde se explica IgV mutated and unmutated (ZAP70): 
# https://watermark.silverchair.com/h81203004944.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAABFwwggRYBgkqhkiG9w0BBwagggRJMIIERQIBADCCBD4GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMvoXw4DTBOiRh7m-MAgEQgIIED90celB-luf785-Ic117_yib0dnmSdfsBbPJySYM8l4wuju-W2TUXbBfG_uFADQnRQx2BYq3PzBIwurXEUfiZBh9OP_EePI0-WLYQftZRYnhbzWmHJH3xqtVvq6NkcSzMn_iXnfsIhZxNZN6eH13cr2S6-D4fLM81aM0w0hI0rPdaKjcHjlsoXI4OCrol5uLiTRTWYTJJCQHL8ZhLqLBbx6ESFRdTT_fECGxOYxh7S2fV7k4KRRK209bPiYRtBAXKOpxxc32sJ23upw77QdZQAaaJgEYgKAuTdQ062B8x6Ibb-CDZXqpn7b23nOsNSNuWpvmjpTUo4C486z0lli25ebQ6MK-sduKchW_sHvNYs8noxFQNCN7Oe67UQBcHBoXaElhfyorRbyFK4zTuwDX1ZizM1-k4P2PyXy1XvlcNyP-lUowLWLpiqd079rajyqh57gJmqvt3GQjKRUuUPIzLkOaGxgQGkiy4G1TONMnE7rdEz2AMyWCOnLKzJo2SOBtooqF51NiMFMUOH2POMOjkjDe3_FZA49Z9WdqPWm8xt2Q3TFrMAK-HQPvvLvqaPhH5HrQSYyiiwRo3Mzq_J6d0tha6JZDeVEB6AgOAWuXO6m4fFugzPp-jr688HtCShiMnFvrPeX5TfLfUfvEzrpaTYJ--y24QM5meoIEwVXDoZx4_Y8fOy88u2cMgDsxKGkKLvFIcrZ4rJE6xXuMh23qEeAfl2bHaJRyty6GLdgjj9YFS89Ni3f8aFe10aSUe2J6Zvyri2ffv1ga06oojjW0_IGkrRHy6UvsZJskQqTSt1_Xgpd1fzm9Pa-Wy-0oF2C4DrEy2OQ2zNRFt537FCk8sPPNCs4e5lo4wP_vRq9YdmOmpW6doogR_YCHBu50R3wXQvo4q_lEqQJfhs1O0OsrXBeoN9PpdII69XpekxioSy5oCGwkch7D_Mnsg2rlqfk2opbKwoV9CU-ZFzQTkfVYG3oN1HcYrDn69whLp0paYkykSKo3iQVRhtUKo5AvPNh7G-Yq8UH_SqOoqKSUWxXYMDRCs_vjEIRr40woOcHp8Ln6er072Q0FxtEt5roR_d7UhcKs-E4JgIEFVTAjCxOqiAEbVMh_16myL9BT0k6RG97yFm0YcD1n5KON0I8jJA2FL1cOtBx0oTVdPQPnTw2noO1qoQTsfabc-4jsjRTygs_ht11B_2yKp86GaTLN5Cat5NygUfCV73ewctLjvdXDWWG-up7AGogQdlq1ZR96ocDPQf-A-1nq4IeProOHgsMjVbxwf7Opt1dQswjh4DryvX6G4G7CtPG5h76OaqTsY9HjhhxaLEzQvZRo6JxmlFdMLCuXo0S03N4J_oyjkSYWMUQ-1x6Bv7Z7_QJqbds3qOw

seurat_annotated = readRDS('Documents/Documents/RT/seurat_annotated.rds')

dim(seurat_annotated)
# 13680   983

Idents(seurat_annotated) = seurat_annotated$orig.ident
VlnPlot(seurat_annotated, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"), ncol = 3, log = F) 

FeaturePlot(seurat_annotated, features = 'nCount_RNA')
FeaturePlot(seurat_annotated, features = 'nFeature_RNA')
FeaturePlot(seurat_annotated, features = 'pct_mt')

Idents(seurat_annotated) = seurat_annotated$annotation_final

seurat_annotated@reductions$umap$UMAP_1

p1 = DimPlot(seurat_annotated, reduction = 'umap', cols = c('darkseagreen', 'palevioletred', 'salmon1'))
p2 = FeaturePlot(seurat_annotated, features = 'ZAP70',order=T, min.cutoff = "q9")
p1 + p2

DotPlot(seurat_annotated, features = c("IL7R", "CCR7", "S100A4", "CD4","CD8A", "CD3E", "CD3G", "GNLY", "NKG7", "MS4A1", "FCGR3A", "MS4A7", "CD14", "LYZ", "FCER1A", "CST3", "PPBP"), col.min = 0) + RotatedAxis()
VlnPlot(seurat_annotated, features = 'ZAP70', log = T, cols = c('darkseagreen', 'palevioletred', 'salmon1'))


seurat_annotated$ZAP_expression <- ifelse(seurat_annotated$RNA$counts['ZAP70', ] >= 1, "ZAP+", "ZAP-")
Idents(seurat_annotated) = seurat_annotated$ZAP_expression


DE_genes_sample <- FindMarkers(object = seurat_annotated, ident.1 = 'ZAP+', ident.2 = 'ZAP-') # Run differential expression analysis
DE_genes_sample = DE_genes_sample[order(-DE_genes_sample[,2], DE_genes_sample[,1] ),]
write.csv(DE_genes_sample, file = "Documents/Documents/RT/DE_genes_63.csv", row.names = T)

DE_genes_sample %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 1.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_genes_sample_pos # aca tengo todos los ZAP70+ genes
write.csv(DE_genes_sample_pos, file = "Documents/Documents/RT/DE_genes_pos63.csv", row.names = T)

DE_genes_sample %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC < -1.5) %>%
  slice_head(n = 500) %>% 
  ungroup() -> DE_genes_sample_neg # aca tengo todos los ZAP70- genes
write.csv(DE_genes_sample_neg, file = "Documents/Documents/RT/DE_genes_neg63.csv", row.names = T)

DE_genes_sample_pos = read.csv('Documents/Documents/RT/DE_genes_pos63.csv')
DE_genes_sample_neg = read.csv('Documents/Documents/RT/DE_genes_neg63.csv')

### GO ENTRICHMENT ANALYSIS
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Entrez IDs
ZAP_entrez_pos <- bitr(DE_genes_sample_pos$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_sample_pos <- enrichGO(gene = ZAP_entrez_pos$ENTREZID,
                                          OrgDb = org.Hs.eg.db,
                                          ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                          pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                          pvalueCutoff = 0.8,
                                          qvalueCutoff = 0.8)
# Visualization
barplot(enrichment_results_sample_pos)
dotplot(enrichment_results_sample_pos)

# Entrez IDs
ZAP_entrez_neg <- bitr(DE_genes_sample_neg$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO analysis
enrichment_results_sample_neg <- enrichGO(gene = ZAP_entrez_neg$ENTREZID,
                                          OrgDb = org.Hs.eg.db,
                                          ont = "BP",  # Biological Process, you can change to "CC" or "MF" if needed
                                          pAdjustMethod = "BH", # Benjamini & Hochberg correction
                                          pvalueCutoff = 0.1,
                                          qvalueCutoff = 0.1)
# Visualization
barplot(enrichment_results_sample_neg)
dotplot(enrichment_results_sample_neg)

###

all_data = readRDS('Documents/Documents/Nadeu2022_NatMed_scRNAseq_data/seurat_objects/5.seurat_list_clustered.rds')
#diznueve = readRDS('Documents/Documents/Nadeu2022_NatMed_scRNAseq_data/seurat_objects/5.seurat_clustered_19.rds')
#doce = readRDS('Documents/Documents/Nadeu2022_NatMed_scRNAseq_data/seurat_objects/5.seurat_clustered_12.rds')
#treseiscinco = readRDS('Documents/Documents/Nadeu2022_NatMed_scRNAseq_data/seurat_objects/5.seurat_clustered_365.rds')
#tresmill = readRDS('Documents/Documents/Nadeu2022_NatMed_scRNAseq_data/seurat_objects/5.seurat_clustered_3299.rds')

#Idents(diznueve) = diznueve$final_clusters
#p1 = DimPlot(diznueve, reduction = 'umap')
#p2 = FeaturePlot(diznueve, features = 'ZAP70',order=T, min.cutoff = "q9")
#p1 + p2

#Idents(doce) = doce$final_clusters
#p1 = DimPlot(doce, reduction = 'umap')
#p2 = FeaturePlot(doce, features = 'ZAP70',order=T, min.cutoff = "q9")
#p1 + p2

#Idents(treseiscinco) = treseiscinco$final_clusters
#p1 = DimPlot(treseiscinco, reduction = 'umap')
#p2 = FeaturePlot(treseiscinco, features = 'ZAP70',order=T, min.cutoff = "q9")
#p1 + p2

#Idents(tresmill) = tresmill$final_clusters
#p1 = DimPlot(tresmill, reduction = 'umap')
#p2 = FeaturePlot(tresmill, features = 'ZAP70',order=T, min.cutoff = "q9")
#p1 + p2

###

annot_diznueve = readRDS('Documents/Documents/RT/Nadeu2022_NatMed_scRNAseq_data/seurat_objects/6.seurat_annotated_19.rds')
dim(annot_diznueve)
# 23326  7284
annot_diznueve$sample <- 19
Idents(annot_diznueve) = annot_diznueve$sample
VlnPlot(annot_diznueve, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"), ncol = 3, log = T) 
saveRDS(annot_diznueve, file = 'Documents/Documents/RT/Nadeu2022_NatMed_scRNAseq_data/seurat_objects/6.seurat_annotated_19.rds')

annot_doce = readRDS('Documents/Documents/RT/Nadeu2022_NatMed_scRNAseq_data/seurat_objects/6.seurat_annotated_12.rds')
dim(annot_doce)
# 23326  5785
annot_doce$sample <- 12
Idents(annot_doce) = annot_doce$sample
VlnPlot(annot_doce, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"), ncol = 3, log = T) 
saveRDS(annot_doce, file = 'Documents/Documents/RT/Nadeu2022_NatMed_scRNAseq_data/seurat_objects/6.seurat_annotated_12.rds')

annot_treseiscinco = readRDS('Documents/Documents/RT/Nadeu2022_NatMed_scRNAseq_data/seurat_objects/6.seurat_annotated_365.rds')
dim(annot_treseiscinco)
# 23326  4685
annot_treseiscinco$sample <- 365
Idents(annot_treseiscinco) = annot_treseiscinco$sample
VlnPlot(annot_treseiscinco, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"), ncol = 3, log = T) 
saveRDS(annot_treseiscinco, file = 'Documents/Documents/RT/Nadeu2022_NatMed_scRNAseq_data/seurat_objects/6.seurat_annotated_365.rds')

annot_tresmill = readRDS('Documents/Documents/RT/Nadeu2022_NatMed_scRNAseq_data/seurat_objects/6.seurat_annotated_3299.rds')
dim(annot_tresmill)
# 23326  6063
annot_tresmill$sample <- 3299
Idents(annot_tresmill) = annot_tresmill$sample
VlnPlot(annot_tresmill, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"), ncol = 3, log = T) 
saveRDS(annot_tresmill, file = 'Documents/Documents/RT/Nadeu2022_NatMed_scRNAseq_data/seurat_objects/6.seurat_annotated_3299.rds')


Idents(annot_diznueve) = annot_diznueve$annotation_final
p1 = DimPlot(annot_diznueve, reduction = 'umap', cols = c('cornflowerblue', 'royalblue', '#F4D35E', 'palevioletred','pink', 'salmon1'))
p2 = FeaturePlot(annot_diznueve, features = 'ZAP70',order=T, min.cutoff = "q9")
p1 + p2
DotPlot(annot_diznueve, features = c("IL7R", "CCR7", "S100A4", "CD4","CD8A", "CD3E", "CD3G", "GNLY", "NKG7", "MS4A1", "FCGR3A", "MS4A7", "CD14", "LYZ", "FCER1A", "CST3", "PPBP"), col.min = 0) + RotatedAxis()

Idents(annot_doce) = annot_doce$annotation_final
p11 = DimPlot(annot_doce, reduction = 'umap', cols = c('cornflowerblue', 'royalblue', 'gold','palevioletred','pink','salmon1', 'violet'))
p22 = FeaturePlot(annot_doce, features = 'ZAP70',order=T, min.cutoff = "q9")
p11 + p22
DotPlot(annot_doce, features = c("IL7R", "CCR7", "S100A4", "CD4","CD8A", "CD3E", "CD3G", "GNLY", "NKG7", "MS4A1", "FCGR3A", "MS4A7", "CD14", "LYZ", "FCER1A", "CST3", "PPBP"), col.min = 0) + RotatedAxis()

Idents(annot_treseiscinco) = annot_treseiscinco$annotation_final
p111 = DimPlot(annot_treseiscinco, reduction = 'umap', cols = c('#F4D35E', 'palevioletred', 'pink', 'violet', 'salmon1'))
p222 = FeaturePlot(annot_treseiscinco, features = 'ZAP70',order=T, min.cutoff = "q9")
p111 + p222
DotPlot(annot_treseiscinco, features = c("IL7R", "CCR7", "S100A4", "CD4","CD8A", "CD3E", "CD3G", "GNLY", "NKG7", "MS4A1", "FCGR3A", "MS4A7", "CD14", "LYZ", "FCER1A", "CST3", "PPBP"), col.min = 0) + RotatedAxis()

Idents(annot_tresmill) = annot_tresmill$annotation_final
p1111 = DimPlot(annot_tresmill, reduction = 'umap', cols = c('darkseagreen','mediumaquamarine','mediumseagreen', 'seagreen','salmon1'))
p2222 = FeaturePlot(annot_tresmill, features = 'ZAP70',order=T, min.cutoff = "q9")
p1111 + p2222
DotPlot(annot_tresmill, features = c("IL7R", "CCR7", "S100A4", "CD4","CD8A", "CD3E", "CD3G", "GNLY", "NKG7", "MS4A1", "FCGR3A", "MS4A7", "CD14", "LYZ", "FCER1A", "CST3", "PPBP"), col.min = 0) + RotatedAxis()

library(gridExtra)
grid.arrange(p1, p2, p11, p22, p111, p222, p1111, p2222, nrow = 4)

markers_19 <- FindAllMarkers(annot_diznueve, only.pos = T)
markers_12 <- FindAllMarkers(annot_doce, only.pos = T)
markers_365 <- FindAllMarkers(annot_treseiscinco, only.pos = T)
markers_3299 <- FindAllMarkers(annot_tresmill, only.pos = T)

markers_19 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)  %>%
  ungroup() -> markers_19
  
markers_12 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)  %>%
  ungroup() -> markers_12

markers_365 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)  %>%
  ungroup() -> markers_365

markers_3299 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)  %>%
  ungroup() -> markers_3299


top_19 <- cut_markers(clusters = levels(markers_19$cluster),markers_19, ntop=100)
top_12 <- cut_markers(clusters = levels(markers_12$cluster),markers_12, ntop=100)
top_365 <- cut_markers(clusters = levels(markers_365$cluster),markers_365, ntop=100)
top_3299 <- cut_markers(clusters = levels(markers_3299$cluster),markers_3299, ntop=100)

## ZAP70 DEG 

annot_diznueve$ZAP_expression <- ifelse(annot_diznueve$RNA$counts['ZAP70', ] >= 1, "ZAP+", "ZAP-")
annot_doce$ZAP_expression <- ifelse(annot_doce$RNA$counts['ZAP70', ] >= 1, "ZAP+", "ZAP-")
annot_treseiscinco$ZAP_expression <- ifelse(annot_treseiscinco$RNA$counts['ZAP70', ] >= 1, "ZAP+", "ZAP-")
annot_tresmill$ZAP_expression <- ifelse(annot_tresmill$RNA$counts['ZAP70', ] >= 1, "ZAP+", "ZAP-")

table(annot_diznueve$ZAP_expression)
# ZAP-  ZAP+ 
# 6743  541 

table(annot_doce$ZAP_expression)
# ZAP-  ZAP+ 
# 5262  523

table(annot_treseiscinco$ZAP_expression)
# ZAP- ZAP+ 
# 4537  148

table(annot_tresmill$ZAP_expression)
# ZAP- ZAP+ 
# 5682  381

Idents(annot_diznueve) = annot_diznueve$ZAP_expression
FeaturePlot(annot_diznueve, features = c("MS4A1","GNLY","CD3E", "CD3G", "CD8A", "IL7R", "CD4", "CD14", "FCGR3A", "LYZ", "PPBP"),order = T) # location in the clusters

markers19 <- FindAllMarkers(annot_diznueve, only.pos = T)
# Save my markers matrix
write.csv(markers19, file = "Documents/Documents/RT/markers19ZAP70.csv", row.names = T)

markers19 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

top <- cut_markers(clusters = levels(markers19$cluster),markers19, ntop=50)

Idents(annot_doce) = annot_doce$ZAP_expression
FeaturePlot(annot_doce, features = c("MS4A1","GNLY","CD3E", "CD3G", "CD8A", "IL7R", "CD4", "CD14", "FCGR3A", "LYZ", "PPBP"),order = T) # location in the clusters

markers12 <- FindAllMarkers(annot_doce, only.pos = T)
# Save my markers matrix
write.csv(markers12, file = "Documents/Documents/RT/markers12ZAP70.csv", row.names = T)

markers12 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

top <- cut_markers(clusters = levels(markers12$cluster),markers12, ntop=50)

Idents(annot_treseiscinco) = annot_treseiscinco$ZAP_expression
FeaturePlot(annot_treseiscinco, features = c("MS4A1","GNLY","CD3E", "CD3G", "CD8A", "IL7R", "CD4", "CD14", "FCGR3A", "LYZ", "PPBP"),order = T) # location in the clusters
markers365 <- FindAllMarkers(annot_treseiscinco, only.pos = T)
# Save my markers matrix
write.csv(markers365, file = "Documents/Documents/RT/markers365ZAP70.csv", row.names = T)

markers365 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

top <- cut_markers(clusters = levels(markers365$cluster),markers365, ntop=50)

Idents(annot_tresmill) = annot_tresmill$ZAP_expression
FeaturePlot(annot_tresmill, features = c("MS4A1","GNLY","CD3E", "CD3G", "CD8A", "IL7R", "CD4", "CD14", "FCGR3A", "LYZ", "PPBP"),order = T) # location in the clusters
markers3299 <- FindAllMarkers(annot_tresmill, only.pos = T)
# Save my markers matrix
write.csv(markers3299, file = "Documents/Documents/RT/markers3299ZAP70.csv", row.names = T)

markers3299%>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

top <- cut_markers(clusters = levels(markers3299$cluster),markers3299, ntop=50)

saveRDS(annot_diznueve, 'Documents/Documents/RT/paciente19.rds')
saveRDS(annot_doce, 'Documents/Documents/RT/paciente12.rds')
saveRDS(annot_treseiscinco, 'Documents/Documents/RT/paciente365.rds')
saveRDS(annot_tresmill, 'Documents/Documents/RT/paciente3299.rds')

###

data_rt <- merge(annot_diznueve, y = c(annot_doce, annot_treseiscinco, annot_tresmill), add.cell.ids = c("S19", "S12", "S365", "S3299"), project = "ALL_SAMPLES_RT")
saveRDS(data_rt, 'RT/data_rt.rds')

data_rt = readRDS('RT/data_rt.rds')

dim(data_rt)
# 23326 23817

data_rt$object <- 1
Idents(data_rt) = data_rt$object
VlnPlot(data_rt, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"), ncol = 3, log = F) 

## 2. NORMALIZATION
Idents(data_rt) = data_rt$
data_rt <- NormalizeData(data_rt, normalization.method = "LogNormalize", scale.factor = 10000)

## 3. IDENTIFY HIGHL VARIABLE FEATURES (feature selection)

data_rt <- FindVariableFeatures(data_rt, selection.method = "vst", nfeatures = 2000) # 15-20% of the number of genes  TRY WITH 800

top10 <- head(VariableFeatures(data_rt), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data_rt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## 4. SCALING - ensures that the magnitudes of gene expression are comparable across all genes, preventing genes with higher magnitudes from dominating the analysis.

data_rt <- ScaleData(data_rt) # This function scales and centers the expression values for each gene across all cells. 

## 5. LINEAR DIMENSIONALITY REDUCTION (PCA)

data_rt <- RunPCA(data_rt, features = VariableFeatures(object = data_rt))

print(data_rt[['pca']], dims = 1:5, nfeatures = 5) # Examine and visualize PRINCIPAL COMPONENT results a few different ways

VizDimLoadings(data_rt, dims = 1:2, reduction = 'pca') # shows the results of the PC's and the genes that are positive and negative. This function visualizes the loadings of features on the first two principal components. It helps you understand which genes contribute the most to the observed variation along these components.

Idents(data_rt) = data_rt$sample
DimPlot(data_rt, reduction = 'pca') + NoLegend() # The position of cells along PC1 and PC2 represents their overall gene expression patterns with respect to the main sources of variation in the dataset. 
### SAVE ###
Idents(data_rt) = data_rt$annotation_final
DimHeatmap(data_rt, dims = 1, cells = 500, balanced = TRUE)
### SAVE ###
DimHeatmap(data_rt, dims = 1:15, cells = 500, balanced = TRUE)

## 6. DETERMINE THE DIMENSIONALITY OF THE dataSET

ElbowPlot(data_rt)

## 7. CLUSTER THE CELLS

data_rt <- FindNeighbors(data_rt, dims = 1:20) #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.

data_rt <- FindClusters(data_rt, resolution = 0.5) #Partition the KNN graph into clusters.

table(data_rt@active.ident)

#  0    1    2    3    4    5    6    7    8    9   10   11   12   13 
# 4221 4118 2873 2595 2396 2169 2014 1762  668  396  255  168  115   67 

head(Idents(data_rt), 5)

## 8. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

data_rt <- RunUMAP(data_rt, dims = 1:20, reduction.name = 'umap.rna')

FeaturePlot(data_rt, features = 'nCount_RNA')
FeaturePlot(data_rt, features = 'nFeature_RNA')
FeaturePlot(data_rt, features = 'pct_mt')

Idents(data_rt) = as.factor(data_rt$sample)
DimPlot(data_rt, reduction = 'umap.rna', label = T, cols = c('cornflowerblue','skyblue4','lightblue','skyblue3'))

Idents(data_rt) = as.factor(data_rt$annotation_final)
DimPlot(data_rt, reduction = 'umap.rna', label = T, cols = c('seagreen','mediumseagreen','darkseagreen','mediumaquamarine','#F4D35E', 'skyblue3', 'royalblue', 'lightblue', 'palevioletred','pink','palevioletred4', 'red3', 'orange', 'sandybrown', 'tan', 'gray'))

Idents(data_rt) = as.factor(data_rt$time_point)
DimPlot(data_rt, reduction = 'umap.rna', label = T, cols = c('seagreen','mediumseagreen','darkseagreen','mediumaquamarine','#F4D35E', 'skyblue3', 'royalblue', 'lightblue', 'palevioletred','pink','palevioletred4', 'red3', 'orange', 'sandybrown', 'tan', 'gray'))

VlnPlot(data_rt, features = 'ZAP70', log = T)

data_rt$ZAP_expression <- ifelse(data_rt$RNA$counts['ZAP70', ] >= 1, "ZAP+", "ZAP-")

# EXPRESSION
FeaturePlot(data_rt, features = 'ZAP70', order  = TRUE, cols = c("grey", "red"), pt.size = 0.5)

# COLOR CELLS THAT ARE ZAP+ 
DimPlot(data_rt, group.by = "ZAP_expression", cols = c("lightgray", "red"), pt.size = 0.1, order = T)

table(data_rt$ZAP_expression)
# ZAP-  ZAP+ 
#  22224  1593

VlnPlot(data_rt, features = 'ZAP70', log = T, c('seagreen','mediumseagreen','darkseagreen','mediumaquamarine','#F4D35E', 'skyblue3', 'royalblue', 'lightblue', 'palevioletred','pink','palevioletred4', 'red3', 'orange', 'sandybrown', 'tan', 'gray'))

## 9. FINDING DIFFERENTIALLY EXPRESSED FEATURES (cluster biomarkers)
Idents(data_rt) = data_rt$annotation_final

markers <- FindAllMarkers(data_rt, only.pos = T)
# Save my markers matrix
write.csv(markers, file = "Documents/Documents/RT/markersZAP70_data_rt.csv", row.names = T)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

DoHeatmap(data_rt, features = top10$gene) + NoLegend()

library(matchSCore2)
top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=100)
top = as.data.frame(top)
write.csv(top, file = "Documents/Documents/RT/top100markers_RT.csv", row.names = T)

top_mydataset = read.csv('Documents/Documents/MULTIOME/top_markers_celltype.csv')
top_mydataset = as.data.frame(top_mydataset)

### SAVE ###
#Most common markers
FeaturePlot(data_rt, features = c("MS4A1","GNLY","CD3E", "CD3G", "CD8A", "IL7R", "CD4", "CD14", "FCGR3A", "LYZ", "PPBP"),order = T) # location in the clusters
VlnPlot(data_rt, features = c("MS4A1","GNLY","CD3E", "CD3G", "CD8A", "CD4", "CD14", "FCGR3A", "LYZ", "PPBP"), slot = "counts", log = TRUE) # RNA counts

DotPlot(data_rt, features = c("GNLY","CD3E", "CD3G", "CD8A", "CD4", "CD14", "FCGR3A", "LYZ", "PPBP")) # RNA counts


## Jaccard Index

## The matchSCore2 function computes the clustering comparison and produce the heatmap table with Jaccard Indexes for each group combination
clustering <- matchSCore2(gene_cl.ref = top ,gene_cl.obs = top_mydataset, ylab = "RT_dataset",xlab = "my_dataset")
## The matchSCore heatmap 
clustering$ggplot

Idents(data_rt) = data_rt$annotation_final
VlnPlot(data_rt, features = 'ZAP70', sort = T)

VlnPlot(data_rt, features = 'ZAP70', sort = T, split.by = 'time_point')

table(data_rt$annotation_final)


