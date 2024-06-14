
## Pipeline for ATAC data
data = readRDS('Documents/Documents/after_meeting/data.NOCD3/final/data_final.rds')

data@active.assay = 'ATAC' # DefaultAssay(data) <- "ATAC"
dim(data)
# 89094 20089

head(rownames(data), 5) # peaks
head(colnames(data), 5) # cells

colnames(data[[]]) # names of the metadata variables
data[["ATAC"]]$counts # the atac counts (peaks)
head(data[[]], 3) # metadata

data
# An object of class Seurat 
# 125695 features across 20089 samples within 2 assays 
# Active assay: ATAC (89094 features, 0 variable features)
# 2 layers present: counts, data
# 1 other assay present: RNA
# 20 dimensional reductions calculated: pca, umap.rna, umap.rna.res.0.5, integrated.cca, umap.cca, integrated.rpca, umap.rpca, harmonyoriginal, umap.harmony.original, integrated.mnn, umap.mnn, harmony, umap.harmony, umap.rna.0.8, umap.rna.1, umap.rna.1.2, umap, umap.no.lymph6, umap.rna.1.no6, umap.harmony.no.lymph6

data_multiome[['ATAC']]

# ChromatinAssay data with 89094 features for 20089 cells
# Variable features: 0 
# Genome: 
# Annotation present: TRUE 
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
# 76    9107   18092   21944   24420  879466 

summary(data$nucleosome_signal)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.5622  0.7086  0.7006  0.8230  2.2747 

summary(data$TSS.enrichment)
#   Min.    1st Qu.  Median  Mean      3rd Qu.  Max. 
# 0.01669  5.15822  5.56083  5.66132  6.00792 52.14785 

#NucleosomeSignal
data$nucleosome_group <- ifelse(data$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
FragmentHistogram(object = data, group.by = 'nucleosome_group')

table(data$nucleosome_group)
#  NS < 2   NS > 2 
#.  20088      2 

#TSS_enrichment
data$high.tss <- ifelse(data$TSS.enrichment > 2.8, 'High', 'Low')
FragmentHistogram(object = data, group.by = 'high.tss')
# TSSPlot(data, group.by = 'high.tss') + NoLegend()

table(data$high.tss)
#  High  Low
# 20045   44 

data <- subset(x = data, subset = nucleosome_signal < 2 & TSS.enrichment > 2.8)

dim(data)
# 89094 20044

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

DimPlot(data, reduction = "lsi", label = F, group.by = 'orig.ident', cols = cols_sam) # latent space, we are visualizing the lsi1/lsi2

data <- RunUMAP(data, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac")
data <- FindNeighbors(object = data, reduction = 'lsi', dims = 2:30)
data <- FindClusters(object = data, verbose = FALSE, algorithm = 3, cluster.name = 'atac_clusters') # SLM algorithm

DimPlot(data, reduction = "umap.atac", group.by = 'atac_clusters',label = T)
DimPlot(data, reduction = "umap.atac", group.by = 'raw_annotation_clusters', label = F, cols = cols_anno)
DimPlot(data, reduction = "umap.atac", group.by = 'orig.ident', label = F, cols = cols_sam)

## Harmony
library(harmony)
data <- RunHarmony(data, "orig.ident", reduction.save = "harmony_atac", assay.use = "ATAC")

data <- RunUMAP(data, dims = 1:20, reduction = 'harmony_atac', reduction.name = 'umap.harmony.atac')

data <- FindNeighbors(data, reduction = "harmony_atac", dims = 1:30, graph.name = "harmony.atac")
data <- FindClusters(data, resolution = 0.5, cluster.name = "harmony_clusters_atac", graph.name = "harmony.atac" )


# basic umaps
DimPlot(data, reduction = 'umap.harmony.atac', group.by = 'new_annotation', label = F, cols = cols_anno)
DimPlot(data, reduction = 'umap.harmony.atac', group.by = 'raw_annotation_clusters', label = F)
DimPlot(data, reduction = 'umap.harmony.atac', group.by = 'prognosis', label = F, cols = cols_group, order = T)
DimPlot(data, reduction = 'umap.harmony.atac', group.by = 'ZAP_expression', label = F, cols = cols_zap, order = T)
DimPlot(data, reduction = 'umap.harmony.atac', group.by = 'orig.ident', label = F, cols = cols_sam)


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

markers_rna <- presto:::wilcoxauc.Seurat(X = data, group_by = 'new_annotation', assay = 'data', seurat_assay = 'RNA') # get RNA markers
markers_motifs <- presto:::wilcoxauc.Seurat(X = data, group_by = 'new_annotation', assay = 'data', seurat_assay = 'chromvar') # get motifs

motif.names <- markers_motifs$feature # extract the motifs names

colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna)) # state the names of the columns in markers_RNA (variables)
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs)) # state the names of the columns in markers_RNA (variables)

markers_rna$gene <- markers_rna$RNA.feature # the last column are the features, now genes from the RNA data
markers_motifs$gene <- ConvertMotifID(data, id = motif.names) # the last column are the motif names, now with the correspoinding ID


# a simple function to implement the procedure above
topTFs <- function(new_annotation, padj.cutoff = 1e-2) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == new_annotation, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
    arrange(-RNA.auc)
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group == new_annotation, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
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
library(monocle)
library(cicero)

## Finding co-accessible networks with Cicero (https://stuartlab.org/signac/articles/cicero)

# Infer potential regulatory interactions between different genomic regions based on single-cell chromatin 
# accessibility data. It helps to identify co-accessible regions, which are genomic loci that tend to have similar 
# accessibility patterns across cells, suggesting potential regulatory relationships. 

# Create the Cicero Object

# Install Cicero
# if (!requireNamespace("remotes", quietly = TRUE))
#   install.packages("remotes")
# remotes::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")

# convert to CellDataSet format and make the cicero object
data.cds <- as.cell_data_set(x = data)
data.cds <- cluster_cells(data.cds)
data.cicero <- make_cicero_cds(data.cds, reduced_coordinates = reducedDims(data.cds)$UMAP.HARMONY.ATAC)

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
library(ggplot2)

clusters <- c("CXCR4loIGHMlo", "CXCR4midIGHMlo", "CXCR4loTXNRD1hi", "CXCR4loIGHMmidMAST4hi", "CXCR4midIGHMloNOTCH2mid", "CXCR4hiIGHMhi|CD27hiMZB1hiXBP1hi")
Idents(data)= as.factor(data$new_annotation)

p = CoveragePlot(object = data, 
                 region = "chr2-97700000-97740000",
                 expression.assay = "RNA",
                 idents = clusters,
                 extend.upstream = 100,
                 extend.downstream = 1000)
p & scale_fill_manual(values = cols_anno)

plot <- CoverageBrowser(object =  data, region = 'ZAP70')

Idents(data)= as.factor(data$prognosis)

clusters <- c("bad.prognosis", "good.prognosis")

clusters <- c("ZAP-", "ZAP+")

p = CoveragePlot(object = data, 
                 region = "ZAP70",
                 expression.assay = "RNA",
                 idents = clusters,
                 extend.upstream = 100,
                 extend.downstream = 1000)

p & scale_fill_manual(values = cols_group)

plot <- CoverageBrowser(object =  data, region = 'ZAP70') 

library(Signac)
library(GenomicRanges)

# obj is an example Seurat object containing a ChromatinAssay

overlap_peaks <- findOverlaps(data, StringToGRanges("chr2-97690000-97750000"))
granges(data)[queryHits(overlap_peaks)]

# chr2 97709623-97710514 

# chr2 97713187-97714102
# chr2 97718305-97719203

# chr2 97740274-97741153

# Find differentially accessible peaks between cell types
DefaultAssay(data) <- 'ATAC'

Idents(data) = as.factor(data$ZAP_expression)
da_peaks_ZAP <- FindAllMarkers(object = data, only.pos = T) # here we calculated the differentially accessible peaks per cluster
head(da_peaks_ZAP)
write.csv(da_peaks_ZAP, file = "Documents/Documents/after_meeting/data.NOCD3/ATAC/da_peaks_ZAP.csv", row.names = T)

da_peaks_ZAP %>%
  group_by(cluster) %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 1.5) %>%
  slice_head(n = 100) %>% 
  ungroup() -> top100_ZAP # aca tengo los 100 mas estadisticamente significativos de cada celltype


Idents(data) = as.factor(data$prognosis)
da_peaks_prognosis <- FindAllMarkers(object = data, only.pos = T) # here we calculated the differentially accessible peaks per cluster
head(da_peaks_prognosis)
write.csv(da_peaks_prognosis, file = "Documents/Documents/after_meeting/data.NOCD3/ATAC/da_peaks_prognosis.csv", row.names = T)

da_peaks_prognosis %>%
  group_by(cluster) %>%
  dplyr::filter(p_val < 0.05) %>%
  slice_head(n = 100) %>% 
  ungroup() -> top100_prognosis # aca tengo los 100 mas estadisticamente significativos de cada celltype

Idents(data) = as.factor(data$new_annotation)
da_peaks_clusters <- FindAllMarkers(object = data, only.pos = T) # here we calculated the differentially accessible peaks per cluster
head(da_peaks_clusters)
write.csv(da_peaks_clusters, file = "Documents/Documents/after_meeting/data.NOCD3/ATAC/da_peaks_clusters.csv", row.names = T)

da_peaks_clusters %>%
  group_by(cluster) %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::filter(avg_log2FC > 1.5) %>%
  slice_head(n = 100) %>% 
  ungroup() -> top100_clusters # aca tengo los 100 mas estadisticamente significativos de cada celltype


library(ggplot2)
library(RColorBrewer)

# Volcano plot of differential peaks for ZAP70
ggplot(da_peaks_ZAP, aes(x = avg_log2FC, y = -log10(p_val), color = cluster)) +
  geom_point(alpha = 0.7) + 
  scale_color_manual(values = cols_zap) +
  geom_point(data = top100_ZAP, aes(x = avg_log2FC, y = -log10(p_val)), color = "black", size = 1) +  # Highlight top 100
  geom_text(data = top100_ZAP[1:2,], aes(label = gene), color = "black", size = 2.5, vjust = -0.8, hjust = -0.0001) + # Add labels for black points
  geom_text(data = top100_ZAP[78:80,], aes(label = gene), color = "black", size = 2.5, vjust = -0.8, hjust = 1) + # Add labels for black points
  labs(x = "Average log2 Fold Change",
       y = "-log10(p-value)",
       title = "Volcano Plot of Differential Accessibility") +
  theme_minimal() +
  facet_wrap(~ cluster, scales = "free")  # Separate plot for each cluster

# Volcano plot of differential peaks for prognosis  

ggplot(da_peaks_prognosis, aes(x = avg_log2FC, y = -log10(p_val), color = cluster)) +
  geom_point(alpha = 0.7) + 
  scale_color_manual(values = cols_group) +
  geom_point(data = top100_prognosis, aes(x = avg_log2FC, y = -log10(p_val)), color = "black", size = 1) +  # Highlight top 100
  geom_text(data = top100_prognosis[101:103,], aes(label = gene), color = "black", size = 2.5, vjust = -0.5, hjust = 0.3) + # Add labels for black points
  geom_text(data = top100_prognosis[1:2,], aes(label = gene), color = "black", size = 2.5, vjust = -0.5, hjust = 0.9) + # Add labels for black points
  labs(x = "Average log2 Fold Change",
       y = "-log10(p-value)",
       title = "Volcano Plot of Differential Accessibility") +
  theme_minimal() +
  facet_wrap(~ cluster, scales = "free")  # Separate plot for each cluster

# Volcano plot of differential peaks for clusters  

ggplot(da_peaks_clusters, aes(x = avg_log2FC, y = -log10(p_val), color = cluster)) +
  geom_point(alpha = 0.7) + 
  scale_color_manual(values = cols_anno) +
  geom_point(data = top100_clusters, aes(x = avg_log2FC, y = -log10(p_val)), color = "black", size = 1) +  # Highlight top 100
  geom_text(data = top100_clusters[101:102,], aes(label = gene), color = "black", size = 2.5, vjust = -0.5) + # Add labels for black points
  geom_text(data = top100_clusters[1:2,], aes(label = gene), color = "black", size = 2.5, vjust = -0.5) + # Add labels for black points
  labs(x = "Average log2 Fold Change",
       y = "-log10(p-value)",
       title = "Volcano Plot of Differential Accessibility") +
  theme_minimal() +
  facet_wrap(~ cluster, scales = "free")  # Separate plot for each cluster

gene = top100_ZAP$gene
# ZAP-
da_peaks_ZAPneg = gene[1:78]
#ZAP+
da_peaks_ZAPpos = gene[79:83]

gene = top100_prognosis$gene
# ZAP-
da_peaks_bad = gene[1:100]
#ZAP+
da_peaks_good = gene[101:200]

gene = top100_clusters$gene

# CXCR4midIGHMloNOTCH2mid
da_peaks_cluster4 = gene[1:100]
#CXCR4hiIGHMhi|CD27hiMZB1hiXBP1hi
da_peaks_cluster5 = gene[101:200]

# The DA peaks do not overlap any region of the ZAP70 gene
find_overlapping_coordinates(da_peaks_bad, "chr2-97000000-97750000")
# character(0)
find_overlapping_coordinates(da_peaks_good, "chr2-97000000-97750000")
# character(0)

find_overlapping_coordinates(da_peaks_ZAPneg, "chr2-97000000-97750000")
# character(0)
find_overlapping_coordinates(da_peaks_ZAPpos, "chr2-97000000-97750000")
# character(0)

find_overlapping_coordinates(da_peaks_cluster4, "chr2-97000000-97750000")
# character(0)
find_overlapping_coordinates(da_peaks_cluster5, "chr2-97000000-97750000")
# character(0)


# We do this for getting the top5 genes of each clusters

# Load necessary packages
library(dplyr)
library(ggplot2)

# Arrange the tibble
top100_clusters_ordered <- top100_clusters %>%
  arrange(p_val_adj, desc(avg_log2FC))

# Subset the first 10 rows
top10_clusters <- top100_clusters_ordered %>%
  slice(1:15)

# Print the subset tibble
print(top10_clusters)

# Create a barplot using ggplot2, ensuring the order is maintained
ggplot(top10_clusters, aes(x = factor(gene, levels = gene), y = p_val_adj)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(x = "Gene", y = "p_val_adj ", title = "Top 15 Peaks by p_val_adj") +
  theme_minimal()

markers_peaks = c('NA', 'GTF2E2', 'NA', 'ZNF827', 'NA', 'NA', 'NA', 'SPATA6L', 'WDFY1', 'NA', 'HSD17B12', 'FOXN2', 'IFNL3', 'HDAC9', 'CHIA')
