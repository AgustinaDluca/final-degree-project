---
title: "RNA_ATAC samples"
author: "Agustina De Luca"
date: "2024-02-19"
output: html_document
---
library(Seurat)
library(SeuratObject)
library(Signac)
library(dplyr)
library(patchwork)
library(PCAtools)
library(devtools)
library(hdf5r)
library(biovizBase)
library(matchSCore2)
library(grDevices)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(Matrix)
library(irlba)

data = readRDS("Documents/Documents/MULTIOME/data_oficial.rds")

## CREATE THE OBJECT WITH THE RNA AND ATAC OMICS!!!!!!!!

### How to load multiome data

### el link para cuando tenga que citar como lo hice (el merge)
## https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/using/aggr?src=pr&lss=none&cnm=&cid=NULL&src=pr&lss=none&cnm=&cid=NULL#requirements

counts <- Read10X_h5("../adeluca/Documents/Documents/MULTIOME/filtered_feature_bc_matrix.h5")
names(counts)

data_init <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA")

# Bash
# tabix -p bed ../adeluca/Documents/Documents/MULTIOME/atac_fragments.tsv.gz

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"
Annotation(data) = annotation

data[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = "Documents/Documents/MULTIOME/atac_fragments.tsv.gz",
  annotation = annotation
)

FeaturePlot(data_init, features = 'nCount_RNA')
FeaturePlot(data_init, features = 'nFeature_RNA')
FeaturePlot(data_init, features = 'percent.mt')

# check for fragment index file
data@assays$ATAC@fragments
# A Fragment object for 25676 cells

data@active.assay = 'RNA'
dim(data)
# 36601 25676
data@active.assay = 'ATAC'
dim(data)
# 89094 25676

saveRDS(data, file = "Documents/Documents/MULTIOME/data_oficial.rds")

## Pipeline for RNA data

data_init@active.assay = 'RNA'
dim(data_init)
# 36601 25676

head(rownames(data), 5) # genes
head(colnames(data), 5) # cells

# Create a column with the samples identifiers
sample_numbers <- sapply(strsplit(rownames(data@meta.data), "-"), function(x) tail(x, 1))
data@meta.data$sample <- sample_numbers

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-") # percentage of mitochondrial genes in each cell

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
#. 0     814    1013    1152    1262    8641

summary(data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   2.869   4.145   4.870   6.124  19.966 

### SELECT THOSE CELLS THAT HAVE A GOOD SCORE REGARDING NUMBER OF FEATURES AND THE PERCENTAGE OF MITOCHONDRIAL CONTAMINATION!!!

data <- subset(data, subset = nFeature_RNA > 800 & percent.mt < 20) # We are going to select the cells we want to be analysed based on the previous results mainly the violin plots from above

dim(data)
# 36601 

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

#.   0    1    2    3    4    5    6    7    8    9 
#. 5763 5733 4928 1324 1248 1088 1020  588  440  244 

head(Idents(data), 5)

## 8. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

data <- RunUMAP(data, dims = 1:20, reduction.name = 'umap.rna')

Idents(data) = as.factor(data$seurat_clusters)
DimPlot(data, reduction = 'umap.rna', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen','cornflowerblue', 'pink', 'royalblue', 'lightblue', 'palevioletred','#FFD700'))

Idents(data) = as.factor(data$sample)
DimPlot(data, reduction = 'umap.rna', label = T, cols = c('darkseagreen','mediumaquamarine','mediumseagreen','seagreen'))

Idents(data) = data$seurat_clusters
VlnPlot(data,features = "nFeature_RNA",log = T, sort = T, c('lightblue','cornflowerblue','palevioletred','royalblue','mediumseagreen','#FFD700','pink', 'seagreen', 'darkseagreen', 'seagreen', 'mediumaquamarine'))

FeaturePlot(data, features = 'nCount_RNA')
FeaturePlot(data, features = 'nFeature_RNA')
FeaturePlot(data, features = 'percent.mt')

## 9. FINDING DIFFERENTIALLY EXPRESSED FEATURES (cluster biomarkers)
Idents(data) = data$seurat_clusters
markers <- FindAllMarkers(data, only.pos = T)

# Save my markers matrix
# write.csv(markers, file = "Documents/Documents/MULTIOME/markers_multiome.csv", row.names = T)
markers = read.csv('Documents/Documents/MULTIOME/markers_multiome.csv')

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

DoHeatmap(data, features = top10$gene) + NoLegend()

library(matchSCore2)
top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=50)
top = as.data.frame(top)
#write.csv(top, file = "Documents/Documents/MULTIOME/top100markers_multiome.csv", row.names = T)
#write.csv(top, file = "Documents/Documents/MULTIOME/top50markers_multiome.csv", row.names = T)

top = read.csv('Documents/Documents/MULTIOME/top100markers_multiome.csv')

### SAVE ###
#Most common markers
FeaturePlot(data, features = c("MS4A1","GNLY","CD3E", "CD3G", "CD8A", "IL7R", "CD4", "CD14", "FCGR3A", "LYZ", "PPBP"),order = T) # location in the clusters
VlnPlot(data, features = c("MS4A1","GNLY","CD3E", "CD3G", "CD8A", "CD4", "CD14", "FCGR3A", "LYZ", "PPBP"), slot = "counts", log = TRUE) # RNA counts
### SAVE ###

# Color each cluster
DimPlot(data, group.by = "seurat_clusters", cols = c('mediumaquamarine','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','mediumseagreen', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'seagreen', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'darkseagreen', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'royalblue', 'grey', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'grey', '#118AB2', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'grey', 'grey', 'cornflowerblue', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'lightblue', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'palevioletred', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey','#F4D35E' ))

top$X0

# "IFNG-AS1"   "RPS26"      "IGKC"       "MAST4"      "JCHAIN"     "TFEC"      
# "HLA-DRB1"   "FUT8"       "IGHD"       "TPD52"      "LINC01374"  "BANK1"     
# "TMSB4X"     "MARCH1"     "GNG7"       "TBC1D9"     "S100A4"     "SAMD12"    
# "RIPOR2"     "CNTNAP2"    "HLA-DRA"    "AL589693.1" "RPS29"      "MT-ND3"    
# "RPL27A"     "MT-ND2"     "MT-ND1"     "RPS27"      "MT-ATP6"    "MT-CO1"    
# "MT-CO3"     "RPL37"      "SLC44A5"    "USP6NL"     "RALGPS2"    "AC108879.1"
# "ALOX5"      "RPL34"      "CAMK4"      "LINC00926"  "RPLP2"      "AL512658.2"
# "DGKG"       "AIM2"       "OSBPL10"    "RPS25"      "TMSB10"     "DIPK1A"    
# "RPS17"      "ALG1L9P"   

top$X1

# "IGLC2"      "IGLC3"      "IGHM"       "GPM6A"      "DLG2"       "ZNF804A"   
# "RBM20"      "BACH1"      "EIF1"       "CXCR4"      "UBC"        "RPL10"     
# "RPS18"      "RPS11"      "EEF1A1"     "TLE4"       "H3F3B"      "RPS3A"     
# "RPS27A"     "RPL13"      "SQSTM1"     "RPS3"       "RPL3"       "TRPS1"     
# "CCDC50"     "HSPA5"      "PTMA"       "RPL9"       "EEF1B2"     "SON"       
# "RPS6"       "HLA-B"      "ZFAND6"     "RPL7"       "RABGAP1L"   "GRIK3"     
# "RPL11"      "BX284613.2" "TRAPPC10"   "TESC"       "TOP1"       "MTMR1"     
# "RPL5"       "ID3"        "UBALD2"     "FOXP1"      "CCNL1"      "TOX2"      
# "RPL13A"     "RBM38"   

top$X2 # NOTCH2 resistance to therapy

# "AC018742.1" "PCDH9"      "TSHZ2"      "SDK2"       "CLNK"       "TRIO"      
# "AC068051.1" "ADK"        "MGAT5"      "TFEC"       "SNED1"      "UFL1-AS1"  
# "SEL1L3"     "FAM126A"    "SNHG14"     "GABPB1-AS1" "ARID5B"     "ITPR2"     
# "SIPA1L1"    "BLK"        "CERS6"      "ADAM28"     "ARHGAP25"   "CLEC2D"    
# "NRCAM"      "PTK2"       "AL390957.1" "LINC00640"  "AL589693.1" "TDRD15"    
# "SKAP1"      "ANKRD44"    "POF1B"      "SMCHD1"     "SMAP2"      "SH3BGRL2"  
# "RGS13"      "IGKC"       "TRAK2"      "KLF7"       "RCSD1"      "LINC01619" 
# "RALGPS2"    "SFMBT1"     "CACNB2"     "RPS12"      "SP100"      "FCRL2"     
# "FNBP4"      "ELOVL5"    

top$X3

# "TEX14"      "AC253572.2" "AC005580.1" "DMD"        "FOSB"       "MT-ND4L"   
# "ENPP2"      "AC007952.4" "AC103591.3" "JUND"       "JUN"        "IQCN"      
# "CACNB2"     "DUSP1"      "SYNE2"      "SGSM1"      "AC090125.1" "LINC00910" 
# "RAB11FIP1"  "ADGRE5"     "AHNAK"      "RGCC"       "KSR2"       "KLF6"      
# "RAPGEF5"    "AC245014.3" "FOS"        "TMEM107"    "DPYD"       "PPP1R9A"   
# "PRCP"       "MSR1"       "AL021155.5" "AC020916.1" "HIP1"       "BCL7A"     
# "EGR1"       "SEPTIN10"   "AC239799.2" "WNT5B"      "Z93241.1"   "TOX"       
# "FGFR1"      "AC012447.1" "IGF2BP3"    "LETM2"      "AL355075.4" "TIMP2"     
# "ANKRD30BL"  "AL356417.3"

top$X4

# "IGLC3"      "IGLC2"      "RBM20"      "GRIK3"      "TRPS1"      "HSPA5"     
# "TESC"       "IRS2"       "DDIT3"      "TOX2"       "U2AF1L5"    "ZC3H12B"   
# "GABRB2"     "APP"        "PTH2R"      "ZNF208"     "ARHGAP42"   "ARSI"      
# "ZNF595"     "VIPR1"      "XBP1"       "MYO3A"      "DLG2"       "ZSCAN18"   
# "CCDC50"     "LRCH1"      "GPM6A"      "BTG3"       "CALR"       "IFI44"     
# "PIM1"       "AC092120.3" "ZNF804A"    "ADNP2"      "ID3"        "MX1"       
# "BACH1"      "UBALD2"     "ZNF793"     "ICOSLG"     "STT3B"      "LINC01725" 
# "CEP68"      "NFKBID"     "UBE2E3"     "PPP1R15A"   "ZNF675"     "PLD5"      
# "AC090125.1" "MTMR1" 

top$X5_0

# "NKG7"       "CCL5"       "PRKCH"      "FYN"        "CD247"      "GNLY"       "AOAH"      
# "STAT4"      "RORA"       "SYTL3"      "GZMH"       "RUNX3"      "ADGRE5"     "TNFAIP3"   
# "IL32"       "PIK3R1"     "BCL11B"     "ITGAL"      "RAP1GAP2"   "NCALD"      "ITGA4"     
# "PRKCQ"      "AAK1"       "ARL4C"      "CST7"       "FYB1"       "GZMM"       "ITGB2"     
# "PRF1"       "GPRIN3"     "SLA2"       "TNIK"       "YES1"       "PDE3B"      "APOL6"     
# "ID2"        "FGFBP2"     "ANXA1"      "ITK"        "KLRK1"      "PPP2R2B"    "CD3D"      
# "PDE4D"      "KLRD1"      "GZMB"       "TRBC1"      "PLCB1"      "IFIT2"      "A2M"       
# "GZMA"       "OASL"       "CD8A"       "SPON2"      "C1orf21"    "NIBAN1"     "LINC01934" 
# "THEMIS"     "SAMD3"      "CD7"        "CCL4"       "CTSW"       "MATK"       "SLAMF7"    
# "CD3E"       "PZP"        "OSBPL5"     "IFIT3"      "ADGRG1"     "AGAP1"      "IFNG"      
# "SH2D2A"     "AC243829.2" "AC015849.1" "GRAP2"      "LINC02384"  "TRGC2"      "CD3G"      
# "KLRF1"      "FCRL6"      "MIAT"       "CCL3"       "DPYD"       "CLIC3"      "SLCO3A1"   
# "S1PR5"      "FCGR3A"     "C12orf75"   "CD2"        "PLAAT4"     "HLA-B"      "ARHGAP10"  
# "ARHGAP26"   "GATA3"      "PRSS23"     "HOPX"       "MYO6"       "IL18RAP"    "COL6A2"    
# "TRERF1"     "IQGAP2"

top$X5_1

# "FYN"        "PRKCH"      "STAT4"      "RORA"       "CCL5"       "SYTL3"      "TC2N"      
# "AOAH"       "CD247"      "PIK3R1"     "DPYD"       "IL7R"       "PDE3B"      "BCL11B"    
# "NKG7"       "TNFAIP3"    "GNLY"       "ITGA4"      "PLCB1"      "IL32"       "AAK1"      
# "MYBL1"      "DUSP2"      "NIBAN1"     "ITK"        "ANK3"       "PRKCQ"      "FYB1"      
# "GZMK"       "TNIK"       "APBA2"      "MGAT4A"     "LINC01934"  "THEMIS"     "INPP4B"    
# "SLCO3A1"    "CST7"       "NCALD"      "CD3D"       "KLRB1"      "CD8A"       "SAMD3"     
# "A2M"        "GRAP2"      "TAFA1"      "PZP"        "PPP2R2B"    "KLRD1"      "RASGRF2"   
# "GZMA"       "SGCD"       "PCAT1"      "CTSW"       "DTHD1"      "FAM169A"    "SYNM"      
# "AGAP1"      "PDE4D"      "TRBC1"      "PAG1"       "ANXA1"      "GZMM"       "GPRIN3"    
# "PTPRM"      "KLRF1"      "ARL4C"      "CFH"        "LINC01871"  "IL12RB2"    "CYTH3"     
# "IL18R1"     "PARP8"      "RUNX3"      "ARHGAP26"   "CD96"       "RNF125"     "CD28"      
# "CD3E"       "MPP7"       "C1orf21"    "PITPNC1"    "TRAT1"      "YES1"       "SORL1"     
# "PRF1"       "PAM"        "SLFN12L"    "TANC2"      "LINC00299"  "NR3C2"      "LRRC8C"    
# "GIMAP5"     "SLA2"       "TNIP3"      "AC010175.1" "TRERF1"     "LINC00612"  "ME1"       
# "LRIG1"      "SPEF2"
top$X6

# "IFNG-AS1"   "MAST4"      "LINC01374"  "JCHAIN"     "LIMCH1"     "ALG1L9P"   
# "SLC44A5"    "MVB12B"     "SAMD12"     "AIM2"       "FUT8"       "CNTNAP2"   
# "HIF1A"      "AL590807.1" "TPD52"      "TFEC"       "CAMK4"      "AC025569.1"
# "FRMD5"      "HSD17B4"    "USP6NL"     "GNG7"       "SSBP2"      "AL512658.2"
# "DGKG"       "AP001825.1" "AC009570.2" "AL365272.1" "PTPRG"      "PDE4DIP"   
# "AC108879.1" "METTL8"     "LAMC1"      "PIP5K1B"    "EPB41L5"    "GLYATL1"   
# "SLC14A1"    "RPS26"      "IGHD"       "BANK1"      "CDK6"       "MAP3K20"   
# "KANSL1L"    "ADAM29"     "DUSP22"     "TMEM156"    "KCNH8"      "DIPK1A"    
# "SERINC5"    "AC116099.1"

top$X7

# "AC018742.1" "PCDH9"      "UFL1-AS1"   "FAM126A"    "CERS6"      "KLF7"      
# "NRCAM"      "PTK2"       "AL390957.1" "LINC00640"  "TDRD15"     "POF1B"     
# "L3MBTL4"    "SH3BGRL2"   "TPRG1"      "RGS13"      "AC068051.1" "LINC01266" 
# "CCDC195"    "UGGT2"      "AC109492.1" "LINC00519"  "CACNB2"     "AC011752.1"
# "SDK2"       "TRAK2"      "SLC9C1"     "DPF3"       "MYO3B"      "ZNF528-AS1"
# "ARHGAP6"    "SLC4A7"     "GAS7"       "MICU3"      "VPS37B"     "KLF3"      
# "MAP4K3-DT"  "TSHZ2"      "KCNMB2-AS1" "ADK"        "HRK"        "TXK"       
# "AC106791.1" "AC009411.2" "AFDN"       "KLHL2"      "EPHA4"      "PCNX2"     
# "RHPN2"      "REEP3"

top$X8_0

# "CD247"      "PRKCH"      "BCL11B"     "PRKCA"      "PDE3B"      "INPP4B"     "ANK3"      
# "FYB1"       "ITK"        "IL32"       "IL7R"       "TNIK"       "GPRIN3"     "PAG1"      
# "NELL2"      "GNAQ"       "CD28"       "CSGALNACT1" "MAL"        "TNFRSF25"   "CD3D"      
# "STAT4"      "TNFAIP3"    "EDA"        "RORA"       "THEMIS"     "AAK1"       "PRKCQ-AS1" 
# "TNFSF8"     "AC139720.1" "BTBD11"     "SORL1"      "TRAT1"      "FYN"        "TRBC1"     
# "DPYD"       "TIAM1"      "GIMAP7"     "PCAT1"      "ATP10A"     "ANXA1"      "FMN1"      
# "RETREG1"    "LINC01891"  "AL136456.1" "NIBAN1"     "CMTM8"      "GATA3"      "LRRC8C"    
# "TC2N"       "ICOS"       "CD2"        "NR3C2"      "CHRM3-AS2"  "ADAM23"     "SLCO3A1"   
# "GRAP2"      "TAFA1"      "TBC1D4"     "PRKCA-AS1"  "AC233976.1" "ANKRD55"    "PDE4D"     
# "LEPROTL1"   "PRKCQ"      "CD3E"       "CD96"       "RASGRF2"    "DHRS3"      "AHR"       
# "ADAMTS17"   "LINC01550"  "MPP7"       "LINC00861"  "RTKN2"      "RNF144A"    "CD4"       
# "TCEA3"      "SYTL3"      "KCNQ1"      "TRERF1"     "DDIT4"      "AC010609.1" "ARL4C"     
# "CD7"        "MYC"        "RNF157"     "CD40LG"     "EIF4E3"     "DNAJB1"     "RCAN3"     
# "DEC1"       "SCML1"      "LRIG1"      "SAMHD1"     "CDC14A"     "HAPLN3"     "DIRC3"     
# "OTULINL"    "MGAT4A"

top$X8_1

# "PRKCH"      "BCL11B"     "CD247"      "GPRIN3"     "STAT4"      "PDE3B"      "SYTL3"     
# "ANK3"       "IFI44"      "PAM"        "CD28"       "OASL"       "FYB1"       "U2AF1L5"   
# "TNFAIP3"    "YES1"       "PRKCQ"      "PDE4D"      "FYN"        "PRKCA"      "GRIK3"     
# "GATA3"      "INPP4B"     "VIPR1"      "ZNF208"     "RORA"       "TRPS1"      "GABRB2"    
# "CCL5"       "PLD5"       "ARSI"       "PIM1"       "APP"        "ARL17B"     "LINC01684" 
# "ZAP70"      "TSPYL2"     "ZC3H12B"    "MYO3A"      "TOX2"       "ZNF331"     "APOL6"     
# "ZNF793"     "ITK"        "BTG3"       "ZNF677"     "ZNF595"     "XBP1"       "TGFBR3"    
# "AOAH"       "PTH2R"      "ARHGAP42"   "LINC01725"  "ZSCAN18"    "IFIT2"      "ARID5A"    
# "NKG7"       "SURF4"      "BCR"        "TRAT1"      "AHDC1"      "SEC24D"     "MED26"     
# "ARL4C"      "IRS2"       "NCALD"      "ZNF880"     "IQGAP2"     "DDIT3"      "USP36"     
# "ADNP2"      "KIFC3"      "PER2"       "ELL2"       "MX1"        "PPFIBP2"    "IFIH1"     
# "NIBAN1"     "UBE2S"      "WARS"       "AC092120.3" "ATG2A"      "AAK1"       "SH3D19"    
# "TICAM1"     "ARHGAP26"   "HERC5"      "DHRS3"      "CRTC2"      "IFIT3"      "ADGRE5"    
# "SLAMF7"     "IFI44L"     "IRF1"       "RBM20"      "TNIK"       "IL2RA"      "PIK3R1"    
# "MYO1D"      "IL7R"

top$X9

# "DPYD"       "AOAH"       "PLXDC2"     "MYO1F"      "NAMPT"      "LRMDA"      "SLC8A1"     "VCAN"       "SLCO3A1"    "FCN1"      
# "TBXAS1"     "CSF3R"      "IRAK3"      "SPI1"       "DOCK5"      "PLCB1"      "DMXL2"      "FAM49A"     "TYMP"       "SLC11A1"   
# "TNFAIP2"    "SAMHD1"     "GNAQ"       "AC020916.1" "MCTP1"      "PLEK"       "BID"        "WDFY3"      "AHR"        "EMILIN2"   
# "HCK"        "PLAUR"      "DENND3"     "IFI30"      "MAML3"      "CTBP2"      "ADGRE2"     "TCF7L2"     "PID1"       "PRAM1"     
# "CLEC7A"     "ARRB1"      "DISC1"      "TTYH3"      "ZMIZ1"      "DUSP6"      "TIMP2"      "CYBB"       "PAG1"       "CPNE8"     
# "S100A9"     "LGALS3"     "FGD4"       "LRP1"       "NLRP3"      "CSF2RA"     "MIR181A1HG" "FOSL2"      "SVIL"       "LYZ"       
# "CD36"       "CLEC12A"    "EPB41L3"    "DAPK1"      "ZNF516"     "RAB20"      "SGK1"       "MICAL2"     "GLT1D1"     "CREB5"     
# "OSCAR"      "CST3"       "TIAM1"      "AC090559.1" "RTN1"       "CD300E"     "KLF4"       "FNIP2"      "LPCAT2"     "LINC00937" 
# "PLXNB2"     "ANPEP"      "KCNQ1"      "RIN2"       "MAFB"       "HK2"        "CYP27A1"    "TREM1"      "IGF2BP2"    "SERPINA1"  
# "STAB1"      "ARHGEF11"   "NIBAN2"     "S100A8"     "IL1B"       "CD300LF"    "LILRB3"     "ST6GALNAC3" "ACTN1"      "CPVL"

# Subsclustering 
Idents(data) <- data$seurat_clusters
data <- FindNeighbors(data, dims = 1:20, graph.name = "graph") #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.

data <- FindSubCluster(data, resolution = 0.1, cluster = '8', graph.name = "graph")

table(data$sub.cluster)
#    0    1    2    3    4    5    6    7   8_0  8_1    9 
#. 5763 5733 4928 1324 1248 1088 1020  588  252  188  244 

Idents(data) <- as.factor(data$sub.cluster)
DimPlot(data, reduction = 'umap.rna', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen','cornflowerblue', 'pink', 'royalblue', 'lightblue', 'palevioletred','violet','#FFD700'))
DimPlot(data, reduction = 'umap.rna', label = F, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen','cornflowerblue', 'pink', 'royalblue', 'lightblue', 'palevioletred','violet','#FFD700'))

markers <- FindAllMarkers(data, only.pos = T)
top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=50)
top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=100)

write.csv(markers, '../adeluca/Documents/Documents/MULTIOME/markers_subclustering.csv')
write.csv(top, '../adeluca/Documents/Documents/MULTIOME/top_markers_subclustering.csv')

markers = read.csv('../adeluca/Documents/Documents/MULTIOME/markers_subclustering.csv')
top = read.csv('../adeluca/Documents/Documents/MULTIOME/top_markers_subclustering.csv')

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%  #Selects the top 10 genes for each cluster.
  ungroup() -> top10

DoHeatmap(data, features = top10$gene) + NoLegend()

data@active.assay = 'RNA'
Idents(data) <- as.factor(data$cell_type)
markers_mydataset <- FindAllMarkers(data, only.pos = T)
top_mydataset <- cut_markers(clusters = levels(markers_mydataset$cluster),markers_mydataset, ntop=100)

write.csv(markers_mydataset, 'Documents/Documents/MULTIOME/markers_celltype.csv')
write.csv(top_mydataset, 'Documents/Documents/MULTIOME/top_markers_celltype.csv')

## Harmony

data <- RunHarmony(data, "sample")
data <- RunUMAP(data, dims = 1:20, reduction = 'pca')
data <- RunUMAP(data, dims = 1:20, reduction = 'harmony', reduction.name = 'umap.harmony')

Idents(data) = data$seurat_clusters
DimPlot(data, reduction = 'umap.harmony', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen','cornflowerblue', 'pink', 'royalblue', 'lightblue', 'palevioletred','#FFD700'))

Idents(data) = as.factor(data$sample)
DimPlot(data, reduction = 'umap.harmony', label = T, cols = c('darkseagreen','mediumaquamarine','mediumseagreen','seagreen'))

Idents(data) = as.factor(data$sub.cluster)
DimPlot(data, reduction = 'umap.harmony', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen','cornflowerblue', 'pink', 'royalblue', 'lightblue', 'palevioletred','violet','#FFD700'))
DimPlot(data, reduction = 'umap.harmony', label = F, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen','cornflowerblue', 'pink','palevioletred1', 'royalblue', 'lightblue', 'palevioletred','violet','#FFD700'))


Idents(data) = as.factor(data$cell_type)
DimPlot(data, reduction = 'umap.harmony', label = T, cols = c('pink', 'palevioletred', 'darkseagreen','mediumaquamarine','mediumseagreen','seagreen','#FFD700', 'salmon1'))

### CELLTYPE ANNOTATION (transfere the data)

## Load the two objects
data = readRDS('Documents/Documents/MULTIOME/data_oficial.rds')
celltypes_data <- readRDS('Documents/Documents/MULTIOME/scRNA   intregrated data/SeuratObjects/DATA.rds')

## Get the rownames
cells_data <- rownames(data@meta.data)
annotation_data <- rownames(celltypes_data@meta.data)

## Remove the S? and add at the end the number of csample it belogs to
cleaned_cell_barcodes <- gsub("^S[1-4]_|-[1]$", "", annotation_data)

cleaned_cell_barcodes[1:5967] <- paste0(cleaned_cell_barcodes[1:5967], "-1")
cleaned_cell_barcodes[5968:13787] <- paste0(cleaned_cell_barcodes[5968:13787], "-2")
cleaned_cell_barcodes[13788:15916] <- paste0(cleaned_cell_barcodes[13788:15916], "-3")
cleaned_cell_barcodes[15917:22247] <- paste0(cleaned_cell_barcodes[15917:22247], "-4")

## Assign row names
rownames(data@meta.data) <- cells_data
rownames(celltypes_data@meta.data) <- cleaned_cell_barcodes

## Data transfere
data$annotation <- celltypes_data$raw_annotation_clusters

Idents(data) = as.factor(data$annotation)
DimPlot(data, reduction = "umap.rna", label = T, cols = c('pink','palevioletred','violet','darkseagreen', 'mediumaquamarine', 'mediumseagreen', 'seagreen', '#F4D35E' ))
DimPlot(data, reduction = "umap.rna", label = F, cols = c('pink','palevioletred','violet','darkseagreen', 'mediumaquamarine', 'mediumseagreen', 'seagreen', '#F4D35E' ))

DimPlot(data, reduction = "umap.harmony", label = T, cols = c('pink','palevioletred','violet','darkseagreen', 'mediumaquamarine', 'mediumseagreen', 'seagreen', '#F4D35E' ))

table(data$annotation, data$sub.cluster)

# 22376 - 22247

### CELLTYPE ANNOTATION (AFTER SEARCHING FOR SPECIFIC MARKERS IN CELLTYPIST)

data <- SetIdent(data, value = as.factor(data$sub.cluster))
cell_type <- c("Leukemic cells 4", "Leukemic cells 2", "Leukemic cells 1", "Leukemic cells 3", "Leukemic cells 2","CD8+ Cytotoxic T cells", "Leukemic cells 4", "Leukemic cells 1", "CD4+ T cells", "Transitioning cells", "Myeloid cells")

cell_type <- cell_type[as.factor(levels(data))]
data$cell_type <- cell_type[as.factor(data$sub.cluster)]

names(cell_type) <- levels(data)

Idents(data) = as.factor(data$cell_type)
DimPlot(data, reduction = "umap.rna", label = T, cols = c('pink','palevioletred','darkseagreen', 'mediumaquamarine', 'mediumseagreen', 'seagreen', '#F4D35E',  'salmon1'))

Idents(data) = as.factor(data$cell_type)
DimPlot(data, reduction = "umap.harmony", label = T, cols = c('pink','palevioletred','darkseagreen', 'mediumaquamarine', 'mediumseagreen', 'seagreen', '#F4D35E','salmon1' ))

# Compare cell_types with annotations

predictions <- table(data$cell_type, data$annotation)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                                                                            low = "gray", high = "darkseagreen4") + ylab("Cell type annotation (RNA)") + xlab("Predicted cell type label (ATAC)") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

####################

## Pipeline for ATAC data
data = readRDS('../adeluca/Documents/Documents/MULTIOME/data_oficial.rds')

data@active.assay = 'ATAC' # DefaultAssay(data) <- "ATAC"
dim(data)
# 89094 22376

head(rownames(data), 5) # peaks
head(colnames(data), 5) # cells

colnames(data[[]])
data[["ATAC"]]$counts
head(data[[]], 3)

data
# An object of class Seurat 
# 125695 features across 22376 samples within 2 assays 
# Active assay: ATAC (89094 features, 0 variable features)
# 2 layers present: counts, data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, umap, harmony

# Remove the "S1_" prefix
new_colnames <- sub("^S1_", "", colnames(data))

# Add the suffix "-1" (or other suffix if needed)
new_colnames <- paste0(new_colnames, "-1") # Change "-1" to the appropriate suffix if needed

# Assign the new column names back to the data
colnames(data) <- new_colnames

# Verify the changes
head(colnames(data))


data[['ATAC']]

# ChromatinAssay data with 89094 features for 22376 cells
# Variable features: 0 
# Genome: 
#   Annotation present: TRUE 
# Motifs present: FALSE 
# Fragment files: 1            --->     data[['ATAC']]@fragments

#We can call granges on a Seurat object with a ChromatinAssay set as the active assay (or on a ChromatinAssay) to see the genomic ranges associated with each feature in the object. 
granges(data)

## 1. QC-filtering of cells - here we will compute some metrics in order to filter out bad quality cells

#a) compute nucleosome signal score per cell [Nucleosomes are basic units of DNA packaging in cells]
data <- NucleosomeSignal(object = data)

# data$nucleosome_percentile --> Percentile rank of the nucleosome signal value.
# data$nucleosome_signal --> The actual value representing the strength or intensity of the nucleosome signal in a genomic region.# Sizes of DNA pieces from each cell
hist(data$nucleosome_signal) # If we see clear peaks in a histogram, it means the DNA is wrapped around proteins (nucleosomes) like it should be. This tells us the data is healthy, the higher nucleosome_signal means better data. 

#b) compute TSS enrichment score per cell
data <- TSSEnrichment(object = data, fast = TRUE) 

# data$TSS.enrichment --> Actual enrichment score around TSS.
# data$TSS.percentile --> Percentile rank of this enrichment score relative to others in the dataset.
hist(data$TSS.enrichment)

# This step might involve enriching the data with information related to transcription start sites (TSS), 
# which are positions in the genome where transcription (the process of making RNA from DNA) starts.

# Given that I do not have metadata from my sample I cannot calculate some quality metrics, such as blacklist ratio and fraction of reads in peaks 
# data$pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters * 100
# data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments

# Plot the relationship between variables stored in the object metadata. 
DensityScatter(data, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

# Plot each metric separatelly ########### 

# Visualize all QC_metrics

VlnPlot(
  object = data,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal', 'nFeature_ATAC'),
  pt.size = 0.1,
  ncol = 5,
)

summary(data$nCount_ATAC)
# Min. 1st Qu.  Median   Mean  3rd Qu.   Max. 
#  30    7133    9082   12950   10513   421291

summary(data$nucleosome_signal)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.5594  0.7051  0.6971  0.8198  2.295 

summary(data$TSS.enrichment)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01547  5.15688  5.55887  5.66263  6.00956 52.14785 

#NucleosomeSignal
data$nucleosome_group <- ifelse(data$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
FragmentHistogram(object = data, group.by = 'nucleosome_group')

table(data$nucleosome_group)
#  NS < 2   NS > 2 
#.  22374      2 

#TSS_enrichment
data$high.tss <- ifelse(data$TSS.enrichment > 2.8, 'High', 'Low')
FragmentHistogram(object = data, group.by = 'high.tss')
# TSSPlot(data, group.by = 'high.tss') + NoLegend()

table(data$high.tss)
#  High  Low
# 22320   56 

data <- subset(x = data, subset = nCount_ATAC > 2000 & nucleosome_signal < 2 & TSS.enrichment > 2.8)

data 

# An object of class Seurat 
# 125695 features across 20855 samples within 2 assays 
# Active assay: ATAC (89094 features, 0 variable features)
# 2 layers present: counts, data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, umap, harmony

VlnPlot(
  object = data,
  features = c('nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal', 'nFeature_ATAC'),
  pt.size = 0.1,
  ncol = 5,
)

### 2. Normalization and linear dimensional reduction (https://stuartlab.org/signac/articles/pbmc_multiomic)
  
# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(data) <- "ATAC"

data <- FindTopFeatures(data, min.cutoff = '5')
data <- RunTFIDF(data)
data <- RunSVD(data)

DepthCor(data) # the first one needs to be removed as it captures sequencing depth (technical variation) rather than biological variation.

# Here we see there is a very strong correlation between the first LSI component and the total number of counts for the cell. 
# We will perform downstream steps without this component as we don’t want to group cells together based on their total sequencing 
# depth, but rather by their patterns of accessibility at cell-type-specific peaks.

## 3. ATAC visualization, as we already have the labels of our cells we do not need to do label tranfere

DimPlot(data, reduction = "lsi", label = T, cols = c('violet','pink','palevioletred','mediumseagreen', 'mediumaquamarine', 'darkseagreen', 'seagreen', '#F4D35E')) # latent space, we are visualizing the lsi1/lsi2

data <- RunUMAP(data, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac")
data <- FindNeighbors(object = data, reduction = 'lsi', dims = 2:30)
data <- FindClusters(object = data, verbose = FALSE, algorithm = 3) # SLM algorithm

Idents(data) = as.factor(data$annotation)
DimPlot(data, reduction = "umap.atac", label = T, cols = c('pink', 'palevioletred', 'violet', 'darkseagreen','mediumaquamarine','mediumseagreen', 'seagreen', '#F4D35E'))

Idents(data) = as.factor(data$cell_type)
DimPlot(data, reduction = "umap.atac", label = T, cols = c('pink', 'palevioletred', 'darkseagreen','mediumaquamarine','mediumseagreen', 'seagreen', '#F4D35E','salmon1'))

Idents(data) = as.factor(data$sample)
DimPlot(data, reduction = "umap.atac", label = T, cols = c('darkseagreen', 'mediumaquamarine', 'mediumseagreen', 'seagreen'))

## Harmony

data <- RunHarmony(data, "sample")
data <- RunUMAP(data, dims = 1:20, reduction = 'lsi')
data <- RunUMAP(data, dims = 1:20, reduction = 'harmony', reduction.name = 'umap.harmony.atac')

Idents(data) = as.factor(data$seurat_clusters)
DimPlot(data, reduction = 'umap.harmony.atac', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen','cornflowerblue', 'pink', 'royalblue', 'lightblue', 'palevioletred','#FFD700'))

Idents(data) = as.factor(data$sub.cluster)
DimPlot(data, reduction = 'umap.harmony.atac', label = T, cols = c('seagreen','mediumaquamarine','darkseagreen','mediumseagreen','cornflowerblue', 'pink', 'royalblue', 'lightblue', 'palevioletred','violet','#FFD700'))

Idents(data) = as.factor(data$sample)
DimPlot(data, reduction = 'umap.harmony.atac', label = T, cols = c('darkseagreen','mediumaquamarine','mediumseagreen','seagreen'))

Idents(data) = as.factor(data$annotation)
DimPlot(data, reduction = 'umap.harmony.atac', label = T, cols = c('pink','palevioletred','violet','darkseagreen','mediumaquamarine','mediumseagreen','seagreen','#F4D35E'))

Idents(data) = as.factor(data$cell_type)
DimPlot(data, reduction = 'umap.harmony.atac', label = T, cols = c('pink','palevioletred','darkseagreen','mediumaquamarine','mediumseagreen','seagreen','#F4D35E','salmon1'))


# Find differentially accessible peaks between cell types

DefaultAssay(data) <- 'ATAC'

Idents(data) = as.factor(data$cell_type)

da_peaks <- FindAllMarkers(object = data, only.pos = T) # here we calculated the differentially accessible peaks per cluster
head(da_peaks)

# Save my da_peaks matrix
write.csv(da_peaks, file = "Documents/Documents/MULTIOME/da_peaks.csv", row.names = T)
da_peaks = read.csv('Documents/Documents/MULTIOME/da_peaks.csv')

da_peaks %>%
  group_by(cluster) %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 1.5) %>%
  slice_head(n = 100) %>% 
  ungroup() -> top100 # aca tengo los 100 mas estadisticamente significativos de cada celltype

write.csv(top100, file = "Documents/Documents/MULTIOME/top100.csv", row.names = T)

library(ggplot2)

# Volcano plot of differential peaks for each cluster
ggplot(da_peaks, aes(x = avg_log2FC, y = -log10(p_val), color = cluster)) +
  geom_point(alpha = 0.7) + 
  scale_color_manual(values = brewer.pal(9, "YlGnBu")) +
  geom_point(data = top100, aes(x = avg_log2FC, y = -log10(p_val)), color = "black", size = 1) +  # Highlight top 50
  labs(x = "Average log2 Fold Change",
       y = "-log10(p-value)",
       title = "Volcano Plot of Differential Accessibility") +
  theme_minimal() +
  facet_wrap(~ cluster, scales = "free")  # Separate plot for each cluster

### top50 DA per cluster
top <- cut_markers(clusters = levels(da_peaks$cluster),da_peaks, ntop=50)
write.csv(top, file = "Documents/Documents/MULTIOME/top50_peaks.csv", row.names = T)

top <- cut_markers(clusters = levels(da_peaks$cluster),da_peaks, ntop=100)
write.csv(top, file = "Documents/Documents/MULTIOME/top100_peaks.csv", row.names = T)

top = read.csv('Documents/Documents/MULTIOME/top100_peaks.csv')

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

## CLUSTERS

lymphoid_CD4 <- ClosestFeature(data, regions = top$`CD4+ T cells`)
# Perform pathway enrichment
enrich_result_lymphoid_CD4 <- enrichGO(gene = lymphoid_CD4$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.06, qvalueCutoff = 0.06)
# Visualize the enriched pathways
dotplot(enrich_result_lymphoid_CD4)
barplot(enrich_result_lymphoid_CD4)

lymphoid_CD8 <- ClosestFeature(data, regions = top$`CD8+ Cytotoxic T cells`)
# Perform pathway enrichment
enrich_result_lymphoid_CD8 <- enrichGO(gene = lymphoid_CD8$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
# Visualize the enriched pathways
dotplot(enrich_result_lymphoid_CD8)
barplot(enrich_result_lymphoid_CD8)

leukemic_1 <- ClosestFeature(data, regions = top$`Leukemic cells 1`)
# Perform pathway enrichment
enrich_result_leukemic_1 <- enrichGO(gene = leukemic_1$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
# Visualize the enriched pathways
dotplot(enrich_result_leukemic_1)
barplot(enrich_result_leukemic_1)

leukemic_2 <- ClosestFeature(data, regions = top$`Leukemic cells 2`)
# Perform pathway enrichment
enrich_result_leukemic_2 <- enrichGO(gene = leukemic_2$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.5, qvalueCutoff = 0.)
# Visualize the enriched pathways
dotplot(enrich_result_leukemic_2)
barplot(enrich_result_leukemic_2)

leukemic_3 <- ClosestFeature(data, regions = top$`Leukemic cells 3`)
# Perform pathway enrichment
enrich_result_leukemic_3 <- enrichGO(gene = leukemic_3$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.06, qvalueCutoff = 0.06)
# Visualize the enriched pathways
dotplot(enrich_result_leukemic_3)
barplot(enrich_result_leukemic_3)

leukemic_4 <- ClosestFeature(data, regions = top$`Leukemic cells 4`)
# Perform pathway enrichment
enrich_result_leukemic_4 <- enrichGO(gene = leukemic_4$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.06, qvalueCutoff = 0.06)
# Visualize the enriched pathways
dotplot(enrich_result_leukemic_4)
barplot(enrich_result_leukemic_4)

myeloid <- ClosestFeature(data, regions = top$`Myeloid cells`)
# Perform pathway enrichment
enrich_result_myeloid <- enrichGO(gene = myeloid$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.06, qvalueCutoff = 0.06)
# Visualize the enriched pathways
dotplot(enrich_result_myeloid)
barplot(enrich_result_myeloid)

transition <- ClosestFeature(data, regions = top$`Transitioning cells`)
# Perform pathway enrichment
enrich_result_transition <- enrichGO(gene = transition$gene_id, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.1, qvalueCutoff = 0.1)
# Visualize the enriched pathways
dotplot(enrich_result_transition)
barplot(enrich_result_transition)

############################

## Examine the accessible regions of each cell to determine enriched motifs. (https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis)

# identify enriched DNA sequence motifs associated with accessible chromatin regions in each cell.

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object

DefaultAssay(data) <- "ATAC"
# JASPAR2020 database
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE)) # fetches matrix data for all matrices in the database matching criteria defined by the named arguments and returns a PFMatrixList object

# Delete peaks that are not in the chromosome sequences from the database. SOLUTION!!
gr <- granges(data)
seq_keep <- seqnames(gr) %in% seqnames(BSgenome.Hsapiens.UCSC.hg38) 
seq_keep <- as.vector(seq_keep)
feat.keep <- GRangesToString(grange = gr[seq_keep])
data[['ATAC']] <- subset(data[['ATAC']], features = feat.keep)
###

motif.matrix <- CreateMotifMatrix(features = granges(data), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
data <- SetAssayData(data, assay = 'ATAC', layer = 'motifs', new.data = motif.object)

seqlevels(granges(data))

# "chr1"       "chr2"       "chr3"       "chr4"       "chr5"       "chr6"       "chr7"      
# "chr8"       "chr9"       "chr10"      "chr11"      "chr12"      "chr13"      "chr14"     
# "chr15"      "chr16"      "chr17"      "chr18"      "chr19"      "chr20"      "chr21"     
# "chr22"      "chrX"       "chrY"       "GL000194.1" "GL000195.1" "GL000205.2" "GL000218.1"
# "GL000219.1" "KI270711.1" "KI270713.1" "KI270721.1" "KI270726.1" "KI270728.1" "KI270731.1"
# "KI270734.1"

# These "KI..." contigs are part of the genome assembly but do not have a specific chromosomal assignment. 
# They are often used in genomic analyses, especially when the exact chromosomal location is not critical 
# or known.

# I have used the EnsDb.Hsapiens.v86 as a reference therefore when working with BSgenome.Hsapiens.UCSC.hg38 it does not find some peaks (seqlevels(BSgenome.Hsapiens.UCSC.hg38))

# calculates a per-cell accessibility score for known motifs
data <- RunChromVAR(object = data, genome = BSgenome.Hsapiens.UCSC.hg38)
# This helps to understanding which transcription factor binding motifs are enriched in accessible chromatin regions across different
# cell types. This helps in elucidating potential regulatory mechanisms controlling gene expression.

# We explore the multimodal dataset to identify key regulators of each cell state.
#devtools::install_github("immunogenomics/presto")
library(presto) # to do a fast differential expression analysis

# - using gene expression data  
# - using chromVAR motif accessibilities

# TFs that return significant results in both tests are identified, and their power as markers of cell types is quantified using the AUC statistic.

# We identify transcription factors that are not only differentially expressed but also have motifs with enriched accessibility

markers_rna <- presto:::wilcoxauc.Seurat(X = data, group_by = 'cell_type', assay = 'data', seurat_assay = 'RNA') # get RNA markers
markers_motifs <- presto:::wilcoxauc.Seurat(X = data, group_by = 'cell_type', assay = 'data', seurat_assay = 'chromvar') # get motifs

motif.names <- markers_motifs$feature # extract the motifs names

colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna)) # state the names of the columns in markers_RNA (variables)
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs)) # state the names of the columns in markers_RNA (variables)

markers_rna$gene <- markers_rna$RNA.feature # the last column are the features, now genes from the RNA data
markers_motifs$gene <- ConvertMotifID(data, id = motif.names) # the last column are the motif names, now with the correspoinding ID


# a simple function to implement the procedure above
topTFs <- function(celltype, padj.cutoff = 1e-2) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
    arrange(-RNA.auc)
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
    arrange(-motif.auc)
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(2, 11, 6, 7)], 
    y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
  )
  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, -avg_auc)
  return(top_tfs)
}

# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis

# identify top markers in NK and visualize
head(topTFs("Leukemic cells 1"), 3)

#          RNA.group  gene   RNA.auc      RNA.pval      motif.group motif.feature motif.auc    motif.pval   avg_auc
# 1 Leukemic cells 1 MEF2C 0.6486485 8.516533e-269 Leukemic cells 1      MA0497.1 0.6250443 1.323206e-171 0.6368464
# 2 Leukemic cells 1 FOXO3 0.5100421  2.026562e-04 Leukemic cells 1      MA0157.2 0.6110338 9.840905e-136 0.5605380
# 3 Leukemic cells 1 GABPA 0.5194673  1.956187e-13 Leukemic cells 1      MA0062.3 0.5999338 2.506584e-110 0.5597005

motif.name <- ConvertMotifID(data, name = 'MEF2C')
gene_plot <- FeaturePlot(data, features = "MEF2C", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(data, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot

head(topTFs("Leukemic cells 2"), 3)

#          RNA.group gene   RNA.auc      RNA.pval      motif.group motif.feature motif.auc motif.pval   avg_auc
# 1 Leukemic cells 2 TCF4 0.5901248 4.275947e-104 Leukemic cells 2      MA0830.2 0.6602243          0 0.6251745
# 2 Leukemic cells 2  REL 0.5407439  1.922248e-28 Leukemic cells 2      MA0101.1 0.6969733          0 0.6188586
# 3 Leukemic cells 2 RELA 0.5085976  5.050146e-07 Leukemic cells 2      MA0107.1 0.7033478          0 0.6059727

head(topTFs("Leukemic cells 3"), 3)

#          RNA.group gene   RNA.auc RNA.pval      motif.group motif.feature motif.auc   motif.pval   avg_auc
# 1 Leukemic cells 3 JUND 0.8297587        0 Leukemic cells 3      MA0491.2 0.6211151 1.312395e-49 0.7254369
# 2 Leukemic cells 3  JUN 0.8033344        0 Leukemic cells 3      MA0488.1 0.6446542 5.447501e-70 0.7239943
# 3 Leukemic cells 3 KLF6 0.8237496        0 Leukemic cells 3      MA1517.1 0.5701416 9.871947e-18 0.6969456

head(topTFs("Leukemic cells 4"), 3)

#          RNA.group  gene   RNA.auc      RNA.pval      motif.group motif.feature motif.auc  motif.pval   avg_auc
# 1 Leukemic cells 4  SPIB 0.5237717  1.573205e-30 Leukemic cells 4      MA0081.2 0.6912288 0.00000e+00 0.6075002
# 2 Leukemic cells 4  SPI1 0.5062044  8.288248e-04 Leukemic cells 4      MA0080.5 0.6804501 0.00000e+00 0.5933273
# 3 Leukemic cells 4 IKZF1 0.6244512 3.899618e-243 Leukemic cells 4      MA1508.1 0.5596707 7.82273e-46 0.5920609

head(topTFs("CD4+ T cells"), 3)

#      RNA.group  gene   RNA.auc      RNA.pval  motif.group motif.feature motif.auc    motif.pval   avg_auc
# 1 CD4+ T cells  ETS1 0.6493504  1.656701e-59 CD4+ T cells      MA0098.3 0.8179099 1.594592e-151 0.7336301
# 2 CD4+ T cells RUNX3 0.7304650 2.996622e-244 CD4+ T cells      MA0684.2 0.7266807  5.398848e-78 0.7285729
# 3 CD4+ T cells  JUNB 0.6493770  4.956992e-60 CD4+ T cells      MA0490.2 0.7449933  8.723831e-91 0.6971851

head(topTFs("CD8+ Cytotoxic T cells"), 3)

#        RNA.group  gene   RNA.auc     RNA.pval            motif.group motif.feature motif.auc    motif.pval
# 1 CD8+ Cytotoxic T cells RUNX3 0.7683140 0.000000e+00 CD8+ Cytotoxic T cells      MA0684.2 0.7891695 8.100721e-162
# 2 CD8+ Cytotoxic T cells  ETS1 0.6123953 5.105194e-44 CD8+ Cytotoxic T cells      MA0098.3 0.8877206 3.092629e-289
# 3 CD8+ Cytotoxic T cells  JUNB 0.5810181 7.173833e-24 CD8+ Cytotoxic T cells      MA0490.2 0.8982632 4.645357e-305

head(topTFs("Transitioning cells"), 3)

#             RNA.group gene   RNA.auc     RNA.pval         motif.group motif.feature motif.auc   motif.pval   avg_auc
# 1 Transitioning cells JUNB 0.7640650 1.208924e-61 Transitioning cells      MA0490.2 0.7281713 3.768948e-27 0.7461182
# 2 Transitioning cells JUND 0.7015581 3.771898e-27 Transitioning cells      MA0491.2 0.7449630 4.861054e-31 0.7232605
# 3 Transitioning cells RBPJ 0.7373181 4.447061e-54 Transitioning cells      MA1116.1 0.6760708 8.256514e-17 0.7066945

head(topTFs("Myeloid cells"), 3)

#       RNA.group  gene   RNA.auc      RNA.pval   motif.group motif.feature motif.auc    motif.pval   avg_auc
# 1 Myeloid cells  SPI1 0.8158413  0.000000e+00 Myeloid cells      MA0080.5 0.9049642 2.719909e-105 0.8604027
# 2 Myeloid cells   FOS 0.7136615 8.197146e-175 Myeloid cells      MA0476.1 0.9647995 4.448177e-138 0.8392305
# 3 Myeloid cells FOSL2 0.6471675  0.000000e+00 Myeloid cells      MA0478.1 0.9614209 4.160922e-136 0.8042942

Tcells_8 = head(topTFs("CD8+ Cytotoxic T cells"), 6)

MotifPlot(object = data,
          motifs = Tcells_8$motif.feature,
          assay = 'ATAC')

Tcells4 = head(topTFs("CD4+ T cells"), 6)

MotifPlot(object = data,
          motifs = Tcells4$motif.feature,
          assay = 'ATAC')

Trans_cells = head(topTFs("Transitioning cells"), 6)

MotifPlot(object = data,
          motifs = Trans_cells$motif.feature,
          assay = 'ATAC')

### ZAP70 multiome script

ZAP.neg = head(topTFs_zap70('ZAP-'),6)

MotifPlot(object = data,
          motifs = ZAP.neg$motif.feature,
          assay = 'ATAC')

ZAP.pos = head(topTFs_zap70('ZAP+'), 6)

MotifPlot(object = data,
          motifs = ZAP.pos$motif.feature,
          assay = 'ATAC')

DefaultAssay(data) <- 'chromvar'

FeaturePlot(object = data,
  features = ZAP.pos$motif.feature,
  min.cutoff = 'q10',
  max.cutoff = 'q90')

# https://stuartlab.org/signac/articles/motif_vignette

library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)

library(cicero)

## Finding co-accessible networks with Cicero (https://stuartlab.org/signac/articles/cicero)

# Infer potential regulatory interactions between different genomic regions based on single-cell chromatin 
# accessibility data. It helps to identify co-accessible regions, which are genomic loci that tend to have similar 
# accessibility patterns across cells, suggesting potential regulatory relationships. 

# Create the Cicero Object

# convert to CellDataSet format and make the cicero object
data.cds <- as.cell_data_set(x = data)
data.cds <- cluster_cells(data.cds)
data.cicero <- make_cicero_cds(data.cds, reduced_coordinates = reducedDims(data.cds)$UMAP.ATAC)

# get the chromosome sizes from the Seurat object
genome <- seqlengths(Annotation(data))

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(data.cicero, genomic_coords = genome.df, sample_num = 100)
head(conns)

# Find cis-co-accessible networks (CCANs)

# Sets of genomic regions that are co-accessible and located on the same chromosome (cis-regulatory interactions). 
# Analyzing CCANs could involve identifying clusters of co-accessible regions within the same chromosome and studying
# their potential regulatory interactions.

ccans <- generate_ccans(conns)
# "Coaccessibility cutoff used: 0.08"
head(ccans)

# Add links to a Seurat object

links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(data) <- links

# We can visualize these links using the CoveragePlot() function, or alternatively we could use the CoverageBrowser() function in an interactive analysis:
# The goal of this step is to associate genomic regions (peaks) with the genes they may be regulating, by analyzing the correlation between gene expression and the accessibility of nearby peaks. This process helps identify potential regulatory regions for each gene.

clusters <- c("Leukemic cells 1", "Leukemic cells 2", "Leukemic cells 3", "Leukemic cells 4", "CD4+ T cells", "CD8+ Cytotoxic T cells", "Myeloid cells", "Transitioning cells")
Idents(data)= as.factor(data$cell_type)
CoveragePlot(object = data, 
             region = "ZAP70",
             features = "ZAP70",
             expression.assay = "RNA",
             idents = clusters,
             extend.upstream = 100,
             extend.downstream = 10000)

plot <- CoverageBrowser(object =  data, region = 'ZAP70')


library(Signac)
library(GenomicRanges)

# obj is an example Seurat object containing a ChromatinAssay

overlap_peaks <- findOverlaps(data, StringToGRanges("chr2-97710000-97750000"))
granges(data)[queryHits(overlap_peaks)]

# chr2 97718305-97719203

# chr2 97740274-97741153
# chr2 97745234-97746123
# chr2 97748298-97749187

## 4. Joint UMAP visualization (RNA & ATAC) 

# build a joint neighbor graph using both assays
data <- FindMultiModalNeighbors(object = data,
                                reduction.list = list("pca", "lsi"), 
                                dims.list = list(1:50, 2:50))

# build a joint UMAP visualization
data <- RunUMAP(object = data,
                nn.name = "weighted.nn",
                reduction.name = "wnn.umap")

Idents(data) = data$annotation
DimPlot(data, label = TRUE, repel = TRUE, reduction = "umap", cols = c('mediumseagreen', 'mediumaquamarine', 'darkseagreen', 'violet','seagreen','pink','#F4D35E','palevioletred')) + NoLegend()

### FIND CLUSTERS
data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
Idents(data) = data$wsnn_res.0.8
DimPlot(data, label = TRUE, repel = TRUE, reduction = "umap", cols = c('mediumseagreen', 'mediumaquamarine', 'darkseagreen', 'violet','seagreen','pink','#F4D35E','palevioletred', 'mediumseagreen', 'mediumaquamarine', 'darkseagreen', 'violet','seagreen','pink','#F4D35E','palevioletred','','pink','#F4D35E','palevioletred')) + NoLegend()

## 5. WNN analysis of 10x Multiome, RNA + ATAC (https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis)

# We can visualize clustering based on gene expression, ATAC-seq, or WNN analysis. The differences are more subtle than in 
# the previous analysis (you can explore the weights, which are more evenly split than in our CITE-seq example), but we find 
# that WNN provides the clearest separation of cell states.

# RNA umap
p1 <- DimPlot(data, reduction = "umap.rna", group.by = "cell_type", label = TRUE, label.size = 2.5, repel = TRUE, cols = c('pink','palevioletred','darkseagreen', 'mediumaquamarine', 'seagreen', 'mediumseagreen', '#F4D35E',  'salmon1'))+ ggtitle("scRNA")
# ATAC umap
p2 <- DimPlot(data, reduction = "umap.atac", group.by = "cell_type", label = TRUE, label.size = 2.5, repel = TRUE, cols = c('pink','palevioletred','darkseagreen', 'mediumaquamarine', 'seagreen', 'mediumseagreen', '#F4D35E',  'salmon1')) + ggtitle("scATAC")
# Weighted Nearest Neighbor umap
p3 <- DimPlot(data, reduction = "wnn.umap", group.by = "cell_type", label = TRUE, label.size = 2.5, repel = TRUE, cols = c('pink','palevioletred','darkseagreen', 'mediumaquamarine', 'seagreen', 'mediumseagreen', '#F4D35E',  'salmon1')) + ggtitle("Weighted Nearest Neighbor")
## plot the three of them together
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

# We can visualize the modality weights that were learned for each cell

## Cluster weights of the RNA assay
VlnPlot(data, features = "RNA.weight", group.by = 'cell_type', sort = TRUE, pt.size = 0.1, cols = c('salmon1','#F4D35E','palevioletred', 'pink', 'seagreen', 'mediumseagreen', 'darkseagreen',  'mediumaquamarine')) + NoLegend()
## Cluster weights of the ATAC assay
VlnPlot(data, features = "ATAC.weight", group.by = 'cell_type', sort = TRUE, pt.size = 0.1, cols = c('mediumaquamarine','darkseagreen','mediumseagreen', 'seagreen', 'pink', 'palevioletred', '#F4D35E',  'salmon1')) + NoLegend()


## SHINY APP
library(ShinyCell)

setwd("Documents/Documents/")
scConf = createConfig(data)
makeShinyApp(data, scConf, gene.mapping = TRUE,
             shiny.title = "Multiome approach")



## IDEAS


## Linking peaks to genes (https://stuartlab.org/signac/articles/pbmc_multiomic)

# For each gene, we can find the set of peaks that may regulate the gene by computing the correlation between gene expression 
# and accessibility at nearby peaks, and correcting for bias due to GC content, overall accessibility, and peak size. 
# Running this step on the whole genome can be time consuming, so here we demonstrate peak-gene links for a subset of genes as an example. 
# The same function can be used to find links for all genes by omitting the genes.use parameter:

DefaultAssay(data) <- "ATAC"

library(BSgenome.Hsapiens.UCSC.hg38)

#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# first compute the GC content for each peak
data <- RegionStats(data, genome = BSgenome.Hsapiens.UCSC.hg38, assay = 'ATAC') # Compute the GC content, region lengths for regions in the assay and add to the feature metadata.

# link peaks to genes
data <- LinkPeaks(object = data, peak.assay = "ATAC", expression.assay = "RNA") # Find peaks that are correlated with the expression of nearby genes.


## Do peak calling

library(Signac)
peaks <- CallPeaks(object = data, group.by = "cell_type", macs2.path = "./Library/Python/3.9/lib/python/site-packages/MACS2/")

CoveragePlot(object = data,
             region = "ZAP70",
             ranges = peaks,
             ranges.title = "MACS2")