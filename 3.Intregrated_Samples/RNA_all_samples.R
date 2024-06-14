---
title: "RNA_all samples"
author: "Agustina De Luca"
date: "2024-02-19"
output: html_document
---
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

data <- readRDS('Documents/Documents/MULTIOME/scRNA   intregrated data/SeuratObjects/DATA.rds')# para el objecto con todas las samples pero no ATAC
saveRDS(data, file = "DATA.rds")

dim(data)
# 36601 22247

head(rownames(data), 5) # genes
head(colnames(data), 5) # cells

data = SetIdent(data, value = data$orig.ident)
data = JoinLayers(data)

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-") # percentage of mitochondrial genes in each cell

# Visualize the different metrics with violin plots
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE) 

VlnPlot(data, features = "percent.mt", ncol = 3, log = F) 

hist(data$percent.mt, main = 'Distribution of mitochondrial genes')
hist(data$nFeature_RNA, main = 'Distribution of the nFeatures')
plot(density(data$nFeature_RNA))

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
# 501    1062    1278    1531    1613    9806 

summary(data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   2.869   4.145   4.870   6.124  19.966 

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

## 6. DETERMINE THE DIMENSIONALITY OF THE DATASET

ElbowPlot(data)

## 7. CLUSTER THE CELLS

data <- FindNeighbors(data, dims = 1:20) #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.

data <- FindClusters(data, resolution = 0.5) #Partition the KNN graph into clusters.

table(data@active.ident)

#    0    1    2    3    4    5    6    7    8    9 
#  6122 4840 4792 1537 1496 1408 1087  646  225   94  

head(Idents(data), 5)

## 8. NON-LINEAR DIMENSINALITY REDUTION (UMAP-tSNE)

data <- RunUMAP(data, dims = 1:20)

Idents(data) = data$seurat_clusters
DimPlot(data, reduction = 'umap', label = T, cols = c('mediumaquamarine','mediumseagreen','seagreen','darkseagreen','palevioletred', '#118AB2', 'cornflowerblue', 'lightblue', '#FFD700', 'royalblue'))

Idents(data) = data$orig.ident
DimPlot(data, reduction = 'umap', label = T, cols = c('mediumseagreen','mediumaquamarine','darkseagreen','seagreen'))

Idents(data) = as.factor(data$cell_type)
DimPlot(data, reduction = 'umap', label = T, cols = c('#F4D35E','violet','pink','pink','#FF6B6B', 'steelblue4', 'lightblue', 'darkseagreen', '#118AB2', '#06D6A0', 'lightblue', 'deepskyblue2', 'seagreen', '#F4D35E', 'palevioletred', 'palevioletred'))

FeaturePlot(data, features = 'nCount_RNA')
FeaturePlot(data, features = 'nFeature_RNA')
FeaturePlot(data, features = 'percent.mt')

## 9. FINDING DIFFERENTIALLY EXPRESSED FEATURES (cluster biomarkers)
markers <- FindAllMarkers(data, only.pos = T)

# Save my markers matrix
#write.csv(markers, file = "Documents/Documents/markers_integrated.csv", row.names = T)
markers = read.csv('Documents/Documents/INTEGRATED DATA/tbls_markers/markers_integrated.csv')

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

DoHeatmap(data, features = top10$gene) + NoLegend()

library(matchSCore2)
top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=100)
top = as.data.frame(top)
#write.csv(top, file = "Documents/Documents/top100markers_integrated.csv", row.names = T)
top = read.csv('Documents/Documents/INTEGRATED DATA/tbls_markers/top100markers_integrated.csv')

### SAVE ###
#Most common markers
FeaturePlot(data, features = c("MS4A1","GNLY","CD3E", "CD3G", "CD8A", "CD4", "CD14", "FCGR3A", "LYZ", "PPBP"),order = T) # location in the clusters
VlnPlot(data, features = c("MS4A1","GNLY","CD3E", "CD3G", "CD8A", "CD4", "CD14", "FCGR3A", "LYZ", "PPBP"), slot = "counts", log = TRUE) # RNA counts
### SAVE ###
VlnPlot(data,features = "nFeature_RNA",log = T, sort = T)

# Color each cluster
DimPlot(data, group.by = "seurat_clusters", cols = c('mediumaquamarine','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','mediumseagreen', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'seagreen', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'darkseagreen', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'palevioletred', 'grey', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'grey', '#118AB2', 'grey', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'grey', 'grey', 'cornflowerblue', 'grey', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'lightblue', 'grey', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', '#F4D35E', 'grey'))
DimPlot(data, group.by = "seurat_clusters", cols = c('grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'royalblue'))

top$X0

# "IGLC2"      "IGLC3"      "GPM6A"      "ZNF804A"    "DLG2"       "RBM20"      "IGHM"       "BACH1"      "TLE4"       "HSPA5"      "ZFAND6"    
# "CCDC50"     "TRPS1"      "SQSTM1"     "SON"        "MTMR1"      "TRAPPC10"   "UBALD2"     "GRIK3"      "UBC"        "EIF1"       "CXCR4"     
# "H3F3B"      "EEF1B2"     "PTMA"       "RPS3A"      "RPL10"      "RPL9"       "RPS18"      "RPS27A"     "RPS11"      "EEF1A1"     "TESC"      
# "RPS3"       "RPS6"       "TOP1"       "RBM38"      "LINC01362"  "RPL3"       "CALR"       "U2AF1L5"    "ID3"        "PPP1R15A"   "TOX2"      
# "PDIA3"      "RPL7"       "BX284613.2" "RABGAP1L"   "RPL11"      "RPL28"      "RPL5"       "CCNL1"      "DDIT3"      "FOXP1"      "HLA-B"     
# "XBP1"       "FAU"        "ACSM3"      "DYRK1A"     "PDCD4"      "RPL13"      "ST3GAL5"    "HNRNPA1"    "JUND"       "PRKCE"      "H1FX"      
# "JUNB"       "MROH9"      "IRS2"       "RPL30"      "MARCH3"     "NPM1"       "UBB"        "PRDM2"      "RPS4X"      "VOPP1"      "PMAIP1"    
# "RSRC2"      "HERPUD1"    "APP"        "GABRB2"     "FMR1"       "GNAS"       "RPL7A"      "FUS"        "PTH2R"      "RPL17"      "MORC3"     
# "RPS5"       "RPS10"      "RPL13A"     "CEP68"      "RPS23"      "GSTP1"      "RPL35A"     "RPS19"      "CREB3L2"    "RPS13"      "RPL21"     
# "ANKRD11"   

top$X1

# "AC018742.1" "PCDH9"      "TSHZ2"      "SDK2"       "AC068051.1" "TFEC"       "UFL1-AS1"   "CLNK"       "TRIO"       "SNED1"      "MGAT5"     
# "ADK"        "FAM126A"    "CERS6"      "SEL1L3"     "NRCAM"      "ARID5B"     "SNHG14"     "LINC00640"  "BLK"        "SIPA1L1"    "ITPR2"     
# "TDRD15"     "AL589693.1" "ADAM28"     "CLEC2D"     "ANKRD44"    "SMCHD1"     "RPS12"      "ARHGAP25"   "SKAP1"      "AL390957.1" "IGKC"      
# "RALGPS2"    "PTK2"       "POF1B"      "SH3BGRL2"   "GABPB1-AS1" "RGS13"      "SMAP2"      "LINC01619"  "SP100"      "SFMBT1"     "RCSD1"     
# "KLF7"       "TRAK2"      "CELF2"      "CUL3"       "FNBP4"      "FCRL2"      "ELOVL5"     "MTSS1"      "MALAT1"     "ZBTB20"     "IMMP2L"    
# "MEF2C"      "SLC4A7"     "SWAP70"     "CACNB2"     "KHDRBS2"    "KLF12"      "PUM2"       "TPRG1"      "FCRL1"      "CCDC91"     "AL592429.2"
# "SESN3"      "EZR"        "ATM"        "GAS7"       "HLA-DQA2"   "MAP3K1"     "FBXL17"     "RB1"        "CDKL1"      "PELI1"      "CHST15"    
# "DNAJC10"    "VPS37B"     "RAPGEF1"    "RPLP1"      "LYN"        "KMT2E"      "KLHL2"      "TMEM123"    "FAM117B"    "RPS14"      "ACTR2"     
# "PDE8A"      "L3MBTL4"    "FCRL5"      "SECISBP2L"  "HERC1"      "LINC01266"  "PRKCB"      "KLF3"       "LINC00519"  "BANK1"      "ARID1B"    
# "FLI1"

top$X2 # NOTCH2 resistance to therapy

# "IFNG-AS1"   "RPS26"      "MAST4"      "IGKC"       "TFEC"       "JCHAIN"     "IGHD"       "FUT8"       "HLA-DRB1"   "GNG7"       "TPD52"     
# "LINC01374"  "S100A4"     "SAMD12"     "TBC1D9"     "MARCH1"     "BANK1"      "AC108879.1" "TMSB4X"     "CNTNAP2"    "DGKG"       "ALOX5"     
# "RIPOR2"     "AL512658.2" "SLC44A5"    "HLA-DRA"    "AL589693.1" "RPL37"      "RPL27A"     "RPS28"      "RPS29"      "RPL37A"     "RPS27"     
# "MT-CO1"     "MT-CO3"     "AIM2"       "USP6NL"     "RPL34"      "RPS21"      "ALG1L9P"    "FCRL1"      "CAMK4"      "RPL38"      "LINC00926" 
# "FAM30A"     "RALGPS2"    "MT-ND3"     "TMSB10"     "RPLP2"      "RPS25"      "FCRL3"      "DIPK1A"     "OSBPL10"    "MT-ATP6"    "RPS17"     
# "CD79A"      "MT-ND1"     "RPL36"      "AL590807.1" "CD52"       "SERINC5"    "HLA-DPA1"   "LIMCH1"     "IGHM"       "RPS8"       "RPL36A"    
# "MALAT1"     "CD74"       "HSD17B4"    "FCRL2"      "STX7"       "MT-ND2"     "RPL30"      "RPL23A"     "RPL39"      "ABCA6"      "RPS12"     
# "MVB12B"     "HIF1A"      "RPLP1"      "RPL31"      "METTL8"     "S100A6"     "VPREB3"     "AUTS2"      "ACTB"       "OAZ1"       "RPL15"     
# "AC009570.2" "NOTCH2"     "RPL27"      "MBNL1"      "DUSP22"     "AC087854.1" "RPS15A"     "LINC01754"  "PDE4DIP"    "LY86"       "RPS16"     
# "KANSL1L"

top$X3

# "TEX14"      "AC253572.2" "AC005580.1" "DMD"        "AC007952.4" "FOSB"       "ENPP2"      "CACNB2"     "IQCN"       "DUSP1"      "MT-ND4L"   
# "PCDH9"      "SGSM1"      "JUN"        "JUND"       "AC103591.3" "LINC00910"  "MT-ATP8"    "AC090125.1" "RGCC"       "SYNE2"      "CD69"      
# "RAPGEF5"    "TMEM107"    "KSR2"       "FOS"        "AHNAK"      "ADGRE5"     "AC245014.3" "RAB11FIP1"  "DPYD"       "PPP1R9A"    "KLF6"      
# "L3MBTL4"    "MSR1"       "PRCP"       "AL021155.5" "CLEC2B"     "SEPTIN10"   "Z93241.1"   "BCL7A"      "HIP1"       "AC020916.1" "AC239799.2"
# "EGR1"       "AC012447.1" "TOX"        "LINGO1"     "IGF2BP3"    "TIMP2"      "FGFR1"      "LETM2"      "ANKRD30BL"  "AL360012.1" "AL356417.3"
# "PTPRK"      "AL139020.1" "OR3A2"      "PLD1"       "ANO5"       "AC020905.1" "MT-ND5"     "MT-ND3"     "MT-ND2"     "MT-ND4"     "MT-CO1"    
# "MT-CO2"     "MT-ATP6"    "KANK2"      "MSI2"       "AKAP12"     "IGF1R"      "MT-ND1"     "MAPK4"      "IGHG3"      "JMJD1C"     "AL355075.4"
# "LPL"        "WNT5B"      "H1FX"       "AC105094.2" "CRY1"       "ATP2A3"     "DTX1"       "PSD3"       "RGS7"       "AP003717.1" "GLDN"      
# "AC008397.1" "PITPNM2"    "WDR74"      "FOXP1"      "AC245060.5" "IGLC1"      "MT-CYB"     "MYO15B"     "ZNF667"     "CLEC4C"     "SPON2"     
# "CROCC"

top$X4

# "PRKCH"     "CD247"     "FYN"       "STAT4"     "RORA"      "CCL5"      "SYTL3"     "NKG7"      "AOAH"      "BCL11B"    "IL32"      "TNFAIP3"  
# "GNLY"      "FYB1"      "DPYD"      "RUNX3"     "AAK1"      "PARP8"     "PDE3B"     "PITPNC1"   "PIK3R1"    "IQGAP2"    "ADGRE5"    "SRGN"     
# "PRKCQ"     "ITK"       "ARL4C"     "TNIK"      "TC2N"      "ITGA4"     "ARHGAP26"  "GZMH"      "ITGAL"     "RAP1GAP2"  "GPRIN3"    "ANXA1"    
# "RASA3"     "SORL1"     "PIK3R5"    "CD96"      "NCALD"     "IL7R"      "CD3D"      "ANK3"      "SLFN12L"   "THEMIS"    "TGFBR3"    "GZMM"     
# "SLCO3A1"   "ITGB2"     "NIBAN1"    "CST7"      "INPP4B"    "MYBL1"     "PDE4D"     "TRBC1"     "PRF1"      "DUSP2"     "LRRC8C"    "CRY1"     
# "SLA2"      "PLCB1"     "ATG2A"     "SPATA13"   "MGAT4A"    "PPP2R2B"   "YES1"      "LINC01934" "ID2"       "CD3E"      "RNF125"    "APOL6"    
# "KLRK1"     "CD8A"      "SAMD3"     "TRERF1"    "GRAP2"     "PAM"       "A2M"       "KLRD1"     "HLA-A"     "PXN"       "PAG1"      "PRKCA"    
# "FGFBP2"    "CD7"       "CD2"       "GZMA"      "TOX"       "MPP7"      "PLAAT4"    "C1orf21"   "CTSW"      "PZP"       "GZMB"      "MIAT"     
# "APBA2"     "ATP10A"    "CD3G"      "MATK"

top$X5

# "IFNG-AS1"   "MAST4"      "JCHAIN"     "TFEC"       "LINC01374"  "FUT8"       "SAMD12"     "USP6NL"     "TPD52"      "CAMK4"      "AIM2"      
# "CNTNAP2"    "LIMCH1"     "AC025569.1" "ALG1L9P"    "MVB12B"     "SLC44A5"    "FRMD5"      "GNG7"       "AL512658.2" "AL590807.1" "RPS26"     
# "KANSL1L"    "HIF1A"      "HSD17B4"    "AC023590.1" "AC009570.2" "PDE4DIP"    "IGHD"       "AL365272.1" "AC108879.1" "DGKG"       "PEG10"     
# "SSBP2"      "GLYATL1"    "NOTCH2"     "AP001825.1" "TBC1D9"     "MAP3K20"    "AL589693.1" "TMEM156"    "SLC14A1"    "METTL8"     "BANK1"     
# "MARCH1"     "FCRL3"      "LAMC1"      "ADAM29"     "FIG4"       "COL18A1"    "KCNH8"      "PIP5K1B"    "PTPRG"      "DIPK1A"     "FCRL1"     
# "DUSP22"     "ENAM"       "BMPR1A"     "EPB41L5"    "SERINC5"    "RHEX"       "AC119396.1" "STX7"       "RAB31"      "LY86"       "ANKUB1"    
# "LINC00278"  "CDK6"       "AC061958.1" "LHFPL2"     "FAM30A"     "DPEP2"      "ALOX5"      "DENND5B"    "AP000787.1" "COL19A1"    "IGKC"      
# "APLP2"      "WDR63"      "HLA-DRB1"   "AP001636.3" "MBNL1"      "WNT3"       "AUTS2"      "AC116099.1" "AC087854.1" "RFTN1"      "PRAG1"     
# "RALGPS2"    "AL392172.2" "SGCE"       "LINC01754"  "AL161781.2" "RIPOR2"     "ABCA6"      "UST"        "MEF2C-AS1"  "PDE5A"      "CR1"       
# "FCRL2"

top$X6

# "IGLC3"      "RBM20"      "GRIK3"      "IGLC2"      "TRPS1"      "HSPA5"      "TESC"       "IRS2"       "U2AF1L5"    "DDIT3"      "ZC3H12B"   
# "TOX2"       "LRCH1"      "XBP1"       "AC092120.3" "ADNP2"      "MX1"        "GABRB2"     "PPFIBP2"    "BTG3"       "APP"        "SEC24D"    
# "ARL17B"     "CALU"       "CASP3"      "ZNF595"     "PIM1"       "RNF144B"    "SLC15A4"    "PTH2R"      "ZC4H2"      "VIPR1"      "ICOSLG"    
# "MYO3A"      "ZNF208"     "Z93930.2"   "ZSCAN18"    "VAMP7"      "ZNF793"     "P2RX1"      "IFI44"      "ARHGAP42"   "VBP1"       "LINC01725" 
# "TCFL5"      "PAM"        "ZNF880"     "LIX1-AS1"   "PLD5"       "ARSI"       "FOXP1-IT1"  "ZNF677"     "CRYM"       "LCN8"       "LINC01684" 
# "DNAJB11"    "GPRASP1"    "IFI44L"     "GARS"       "CALR"       "DERL2"      "ZNF675"     "IL2RA"      "PIM2"       "ZNF37A"     "PHEX"      
# "ESR2"       "APBB2"      "SURF4"      "ATP2A2"     "TNFRSF13B"  "JADE3"      "HERC5"      "ZBTB21"     "CMTM6"      "AP002370.2" "GAS6-AS1"  
# "LINC02133"  "PFKL"       "SH3D19"     "IL6ST"      "AL138828.1" "RSL1D1"     "AL109930.1" "ELL"        "ZNF331"     "CCDC50"     "CRAMP1"    
# "LANCL2"     "TSPYL2"     "GRIK1"      "IGHV3-53"   "ACSM1"      "ST3GAL5"    "XRRA1"      "MAML1"      "MYO1D"      "RRBP1"      "TIMM44"    
# "OSTN-AS1"  

top$X7

# "AC018742.1" "PCDH9"      "UFL1-AS1"   "CERS6"      "NRCAM"      "PTK2"       "AL390957.1" "TDRD15"     "LINC00640"  "POF1B"      "SH3BGRL2"  
# "TPRG1"      "RGS13"      "FAM126A"    "KLF7"       "AC068051.1" "L3MBTL4"    "AC011752.1" "LINC01266"  "SDK2"       "UGGT2"      "CCDC195"   
# "SLC9C1"     "LINC00519"  "TRAK2"      "GAS7"       "CACNB2"     "ADK"        "TSHZ2"      "DPF3"       "AC109492.1" "MICU3"      "SLC4A7"    
# "MYO3B"      "ARHGAP6"    "AC009411.2" "SNED1"      "VPS37B"     "HRK"        "KLF3"       "AFDN"       "KCNMB2-AS1" "AL358332.1" "CLNK"      
# "EPHA4"      "RHPN2"      "SIPA1L1"    "TRIO"       "HHAT"       "MGAT5"      "KLHL2"      "PCNX2"      "AC106791.1" "PDE8A"      "MAP4K3-DT" 
# "TXK"        "JADE1"      "CHST15"     "AL589693.1" "ZNF528-AS1" "SEL1L3"     "NRIP1"      "ITPR2"      "PTPRO"      "SNHG14"     "AC009522.1"
# "REEP3"      "ARID5B"     "AL109930.1" "GABPB1-AS1" "U2AF1"      "LYPLAL1"    "GPHN"       "DNAJC10"    "FCRL2"      "ITGAX"      "AL591686.2"
# "TFEC"       "AC025569.1" "C3orf35"    "FBXO10"     "AC002377.1" "ZC2HC1A"    "AC079447.1" "SMCHD1"     "SLC39A10"   "ANKRD44"    "SMAP2"     
# "TCP11L2"    "FCRL5"      "FOCAD"      "TMEM241"    "ARHGAP22"   "MIR646HG"   "SFMBT2"     "LINC01619"  "ADAM28"     "TTC28"      "HOMER3"    
# "HIP1R"

top$X8

# "DPYD"       "PLXDC2"     "AOAH"       "LRMDA"      "SLC8A1"     "VCAN"       "CSF3R"      "FCN1"       "SLCO3A1"    "AC020916.1" "TYMP"      
# "FAM49A"     "TBXAS1"     "DOCK5"      "IRAK3"      "DMXL2"      "PLCB1"      "SLC11A1"    "GNAQ"       "BID"        "TNFAIP2"    "PLEK"      
# "GAS7"       "MCTP1"      "PLAUR"      "AHR"        "HCK"        "EMILIN2"    "IFI30"      "ADGRE2"     "WDFY3"      "CLEC7A"     "TCF7L2"    
# "DENND3"     "MAML3"      "PRAM1"      "TTYH3"      "NLRP3"      "PID1"       "CTBP2"      "ARRB1"      "CSF2RA"     "DISC1"      "DUSP6"     
# "CYBB"       "ZMIZ1"      "FOSL2"      "S100A9"     "LRP1"       "FGD4"       "OSCAR"      "CREB5"      "CPNE8"      "LGALS3"     "ICAM1"     
# "MIR181A1HG" "SGK1"       "GLT1D1"     "RAB20"      "EPB41L3"    "CD36"       "LYZ"        "NINJ1"      "MICAL2"     "CD300E"     "SVIL"      
# "DAPK1"      "HK2"        "ZNF516"     "KLF4"       "IL1B"       "CLEC12A"    "CYP27A1"    "IGF2BP2"    "FNIP2"      "RTN1"       "ANPEP"     
# "SERPINA1"   "TREM1"      "STAB1"      "LINC00937"  "FCAR"       "NECTIN2"    "LILRB3"     "KCNQ1"      "AC090559.1" "PLXNB2"     "CST3"      
# "MAFB"       "ZNF385A"    "RIN2"       "SIRPA"      "INSR"       "S100A8"     "LPCAT2"     "CD300LB"    "CPVL"       "APOBEC3A"   "LINC00877" 
# "EEPD1"     

top$X9

# "ENPP2"      "AC005580.1" "EGR1"       "PPP1R9A"    "TIMP2"      "MSR1"       "PTPRD"      "BCL7A"      "CTNNA3"     "LETM2"      "LRP1B"     
# "IGF2BP3"    "IGF1R"      "RGS7"       "SEPTIN10"   "ANKRD30BL"  "MIR3681HG"  "RBMS3"      "CADM2"      "GRID2"      "NAALADL2"   "ERBB4"     
# "ROBO2"      "PSD3"       "PTPRK"      "RBFOX1"     "DCC"        "AC007402.1" "ZNF385D"    "DPP10"      "CDH18"      "OR3A2"      "SGCZ"      
# "AC060765.2" "PCDH15"     "CSMD3"      "MAGI2"      "CNTN5"      "GALNTL6"    "AC113414.1" "PRKG1"      "GRM7"       "AL391117.1" "LRRC4C"    
# "AC011287.1" "ANKS1B"     "CDH12"      "TRPM3"      "CNTN4"      "LSAMP"      "AC116362.1" "NRG3"       "GPC6"       "PARD3B"     "SOX5"      
# "ANO5"       "AC013287.1" "KCNIP4"     "AC109466.1" "DGKB"       "NKAIN2"     "AL589740.1" "NRG1"       "TENM3"      "AC016766.1" "CNBD1"     
# "MGAT4C"     "ADGRV1"     "LAMA2"      "CCSER1"     "LINC02055"  "IL1RAPL1"   "NRXN1"      "ZFPM2"      "AC078923.1" "AC034268.2" "FRMPD4"    
# "SNTG1"      "AC034114.2" "AC068633.1" "LINC01322"  "NRXN3"      "KIAA1217"   "PARD3"      "ROBO1"      "CECR2"      "GTF2IRD1"   "OBI1-AS1"  
# "AC011246.1" "CASC15"     "DLC1"       "AL157944.1" "EPHA6"      "NELL1"      "SYT1"       "CTNNA2"     "TENM2"      "AC007389.1" "CCDC26"    
# "RIMS2"

# Subsclustering 
Idents(data) <- data$seurat_clusters
data <- FindNeighbors(data, dims = 1:20, graph.name = "graph") #Euclidean distance in the PCA space to define edges between cells with similar feature expression patterns.
data <- FindSubCluster(data, resolution = 0.2, cluster = '4', graph.name = "graph")
table(data$sub.cluster)

Idents(data) <- as.factor(data$sub.cluster)

DimPlot(data, reduction = 'umap', label = T, cols = c('mediumaquamarine','mediumseagreen','seagreen','darkseagreen','palevioletred','pink', 'violet', '#118AB2', 'cornflowerblue', 'lightblue', '#F4D35E', 'royalblue'))
DimPlot(data, reduction = 'umap', label = F, cols = c('mediumaquamarine','mediumseagreen','seagreen','darkseagreen','palevioletred','pink', 'violet', '#118AB2', 'cornflowerblue', 'lightblue', '#F4D35E', 'royalblue'))

markers <- FindAllMarkers(data, only.pos = T)
top <- cut_markers(clusters = levels(markers$cluster),markers, ntop=50)

#write.csv(top1, 'top_markers_subclustering.csv')
markers = read.csv('Documents/Documents/INTEGRATED DATA/tbls_markers/markers_subclustering_data.csv')
top = read.csv('Documents/Documents/INTEGRATED DATA/tbls_markers/top_markers_subclustering_data.csv')

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>% #Selects the top 10 genes for each cluster.
  ungroup() -> top10

DoHeatmap(data, features = top10$gene) + NoLegend()

## Harmony

data <- RunHarmony(data, "orig.ident")
data <- RunUMAP(data, dims = 1:20, reduction = 'harmony')

data_harmony <- readRDS('Documents/Documents/DATA_harmony.rds')

Idents(data_harmony) = data_harmony$seurat_clusters
DimPlot(data_harmony, reduction = 'umap', label = T, cols = c('mediumaquamarine','mediumseagreen','seagreen','darkseagreen','palevioletred', '#118AB2', 'cornflowerblue', 'lightblue', '#FFD700', 'royalblue', 'green', 'yellow'))

Idents(data_harmony) = data_harmony$orig.ident
DimPlot(data_harmony, reduction = 'umap', label = T, cols = c('mediumseagreen','mediumaquamarine','darkseagreen','seagreen'))

Idents(data_harmony) = as.factor(data_harmony$cell_type)
DimPlot(data_harmony, reduction = 'umap', label = T, cols = c('#F4D35E','violet','pink','pink','#FF6B6B', 'steelblue4', 'lightblue', 'darkseagreen', '#118AB2', '#06D6A0', 'lightblue', 'deepskyblue2', 'seagreen', '#F4D35E', 'palevioletred', 'palevioletred'))

Idents(data_harmony) = as.factor(data_harmony$sub.cluster)
DimPlot(data_harmony, reduction = 'umap', label = T, cols = c('mediumaquamarine','mediumseagreen','seagreen','darkseagreen','palevioletred','pink', 'violet', '#118AB2', 'cornflowerblue', 'lightblue', '#F4D35E', 'royalblue'))

### CELLTYPE ANNOTATION (AFTER SEARCHING FOR SPECIFIC MARKERS IN CELLTYPIST)

data <- SetIdent(data, value = as.factor(data$sub.cluster))
raw_annotation_clusters <- c("Leukemic cells 2", "Leukemic cells 1", "Leukemic cells 4", "Leukemic cells 3", "CD8+ Cytotoxic T cells","CD4+ T cells", "Cytotoxic NK", "Leukemic cells 4", "Leukemic cells 2", "Leukemic cells 1", "Myeloid cells", "Leukemic cells 3")

raw_annotation_clusters <- raw_annotation_clusters[as.factor(levels(data))]
data$raw_annotation_clusters <- raw_annotation_clusters[as.factor(data$sub.cluster)]

names(raw_annotation_clusters) <- levels(data)

data <- SetIdent(data, value = as.factor(data$sub.cluster))

DimPlot(data, reduction = "umap", label = T, pt.size = 0.5,  cols = c('violet','pink','palevioletred','mediumseagreen', 'mediumaquamarine', 'darkseagreen', 'seagreen', '#F4D35E')) + NoLegend()
DimPlot(data, reduction = "umap", label = F, pt.size = 0.5,  cols = c('violet','pink','palevioletred','mediumseagreen', 'mediumaquamarine', 'darkseagreen', 'seagreen', '#F4D35E')) + NoLegend()

## SHINY APP
library(ShinyCell)

setwd("/home/adeluca/Documents/")
scConf = createConfig(data)
makeShinyApp(data, scConf, gene.mapping = TRUE,
             shiny.title = "Samples integrated")
