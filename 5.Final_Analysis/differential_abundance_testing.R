####################
## Differential abundance testing - sccomp_glm 

# (https://www.bioconductor.org/packages/release/bioc/vignettes/sccomp/inst/doc/introduction.html#installation)
# (https://www.bioconductor.org/packages/release/bioc/manuals/sccomp/man/sccomp.pdf)

#if (!requireNamespace("BiocManager")) install.packages("BiocManager")
#BiocManager::install("sccomp")

library(sccomp)

data$prognosis <- ifelse(data$orig.ident %in% c('S2', 'S3'), "bad.prognosis", "good.prognosis")

data$prognosis = as.factor(data$prognosis)
data$raw_annotation_clusters = as.factor(data$raw_annotation_clusters)
data$orig.ident = as.factor(data$orig.ident)

table(data$prognosis)
# bad.prognosis good.prognosis 
#.   8462          11627 

table(data$orig.ident)
# S1   S2   S3   S4 
#.5455 6947 1515 6172 

library(sccomp)
res = data |> 
  sccomp_glm( formula_composition = ~ prognosis, # generate a column in my dataset separating the cells into good|bad diagnosis depending on the ZAP70 expression and clinical data
              formula_variability = ~ 1, 
              percent_false_positive = 5, 
              .sample = orig.ident, # each samples
              .cell_group = raw_annotation_clusters) # my final clusters annotated

plots = plot_summary(res) 

plots$boxplot[[1]]

plots$boxplot
ggsave("after_meeting/data.NOCD3/boxplots.pdf", width = 15, height = 8)

plots$credible_intervals_2D
ggsave("after_meeting/data.NOCD3/credible_intervals_2D.pdf", width = 15, height = 8)

plots$credible_intervals_1D
ggsave("after_meeting/data.NOCD3/credible_intervals_1D.pdf", width = 15, height = 8)

## Now we will try the same but with the MILO tool (another one to do differential abundance)

# https://rawcdn.githack.com/MarioniLab/miloR/7c7f906b94a73e62e36e095ddb3e3567b414144e/vignettes/milo_gastrulation.html#5_Finding_markers_of_DA_populations (pipeline)
# https://github.com/MarioniLab/miloR

# BiocManager::install(c("miloR", "scater", "scran"))
library(miloR)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(scater)
library(scran)
library(SummarizedExperiment)

data_sce <- as.SingleCellExperiment(data)
data_milo <- Milo(data_sce)
dim(data_milo)
#dim: 36601 20044 

#data_milo <- buildGraph(data_milo, k = 30, d = 30, reduced.dim = "pca.corrected")
miloR::graph(data_milo) <- miloR::graph(buildFromAdjacency(agus, k=20))

#data_milo <- makeNhoods(data_milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "pca.corrected")

dec_data <- modelGeneVar(data_sce)
fit_data <- metadata(dec_data)
plot(fit_data$mean, fit_data$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression") # A mi este plot me da como una "diagonal" no muy pendiente.
hvgs_data <- getTopHVGs(dec_data, n=3000)

data_milo <- makeNhoods(data_milo, k=30, d=50, prop = 0.05, refined=TRUE)
plotNhoodSizeHist(data_milo)

data_milo <- countCells(data_milo, meta.data = as.data.frame(colData(data_milo)), sample="orig.ident")
head(nhoodCounts(data_milo))

data_design <- data.frame(colData(data_milo))[,c("orig.ident", "prognosis")]

## Convert batch info from integer to factor
data_design <- distinct(data_design)
rownames(data_design) <- data_design$orig.ident
data_design$sample <- NULL
data_design
data_design$prognosis <- factor(data_design$prognosis)

data_milo <- calcNhoodDistance(data_milo, d=50, reduced.dim = "PCA")
data_milo@.k <- 30

data_res <- testNhoods(data_milo, design = ~ prognosis, design.df = data_design, reduced.dim="UMAP.HARMONY.NO.LYMPH6", fdr.weighting="graph-overlap")

head(data_res)

data_res %>%
  arrange(SpatialFDR) %>%
  head() 

ggplot(data_res, aes(PValue)) + geom_histogram(bins=50)

ggplot(data_res, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

data_milo <- buildNhoodGraph(data_milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(data_milo, dimred = "UMAP.HARMONY.NO.LYMPH6", colour_by = "prognosis", text_by = "new_annotation", text_size = 3) +
  guides(fill = "none") +
  scale_colour_manual(values = cols_group)

nh_graph_pl <- plotNhoodGraphDA(data_milo, data_res, layout="UMAP.HARMONY.NO.LYMPH6", alpha=0.05) 

umap_pl + nh_graph_pl +  plot_layout(guides="collect")

data_res <- annotateNhoods(data_milo, data_res, coldata_col = "new_annotation")
head(data_res)

ggplot(data_res, aes(new_annotation_fraction)) + geom_histogram(bins=50)

data_res$new_annotation_fraction <- ifelse(data_res$new_annotation_fraction < 0.8, "Mixed", data_res$new_annotation)

data_res$new_annotation
plotDAbeeswarm(data_res, group.by = "new_annotation") # DA results per cluster

## Automatic grouping of neighbourhoods

data_milo <- buildNhoodGraph(data_milo)

plotDAbeeswarm(groupNhoods(data_milo, data_res, max.lfc.delta = 0.5), group.by = "NhoodGroup") + ggtitle("max LFC delta=0.5")
plotDAbeeswarm(groupNhoods(data_milo, data_res, max.lfc.delta = 1) , group.by = "NhoodGroup") + ggtitle("max LFC delta=1")
plotDAbeeswarm(groupNhoods(data_milo, data_res, max.lfc.delta = 2)   , group.by = "NhoodGroup") + ggtitle("max LFC delta=2")
plotDAbeeswarm(groupNhoods(data_milo, data_res, max.lfc.delta = 3)   , group.by = "NhoodGroup") + ggtitle("max LFC delta=3")
plotDAbeeswarm(groupNhoods(data_milo, data_res, max.lfc.delta = 4)   , group.by = "NhoodGroup") + ggtitle("max LFC delta=3")

## Final values

plotDAbeeswarm(groupNhoods(data_milo, data_res, max.lfc.delta = 4, overlap=1), group.by = "NhoodGroup") + ggtitle("max.lfc.delta=3 | overlap=1")
plotDAbeeswarm(groupNhoods(data_milo, data_res, max.lfc.delta = 3, overlap=1), group.by = "NhoodGroup") + ggtitle("max.lfc.delta=3 | overlap=1")


## PLOTS

data_res <- groupNhoods(data_milo, data_res, max.lfc.delta = 8, overlap = 1)

one = plotNhoodGroups(data_milo, data_res, layout="UMAP.HARMONY.NO.LYMPH6")
two = plotDAbeeswarm(data_res, group.by = "NhoodGroup")

one + two +  plot_layout(guides="collect")

## plot clusters/groups
plotDAbeeswarm(data_res, group.by = 'NhoodGroup') +
  facet_grid(new_annotation~., scales="free", space="free")

## Finding gene signatures for neighbourhoods
## Exclude zero counts genes
keep.rows <- rowSums(logcounts(data_milo)) != 0
data_milo <- data_milo[keep.rows, ]

## Find HVGs
dec <- modelGeneVar(data_milo)
hvgs <- getTopHVGs(dec, n=2000)
head(hvgs)

nhood_markers <- findNhoodGroupMarkers(data_milo, data_res, subset.row = hvgs, 
                                       aggregate.samples = TRUE, sample_col = "orig.ident")

head(nhood_markers)

## Find markers from different Nhoods

nhood_markers2 <- findNhoodGroupMarkers(data_milo, data_res, subset.row = hvgs, 
                                        aggregate.samples = TRUE, sample_col = "orig.ident",
                                        subset.groups = c("2"))
head(nhood_markers2)

#             logFC_2   adj.P.Val_2   GeneID
#AL359232.1 0.03048807   0.9987429 AL359232.1
#AMN1       0.06107171   0.9987429       AMN1
#TTC37      0.08767465   0.9987429      TTC37
#AC078785.1 0.02113458   0.9987429 AC078785.1
#AC112206.2 0.01688290   0.9987429 AC112206.2
#ATRN       0.05026107   0.9987429       ATRN

nhood_markers4 <- findNhoodGroupMarkers(data_milo, data_res, subset.row = hvgs, 
                                        aggregate.samples = TRUE, sample_col = "orig.ident",
                                        subset.groups = c("4"))
head(nhood_markers4)

#            logFC_4    adj.P.Val_4   GeneID
#MT-ND3    -0.08525609   0.9997964    MT-ND3
#DTWD2      0.04122808   0.9997964     DTWD2
#MT-CO2    -0.08090530   0.9997964    MT-CO2
#ATR        0.07809195   0.9997964       ATR
#MCPH1-AS1  0.04213670   0.9997964 MCPH1-AS1
#MYO3A     -0.04687951   0.9997964     MYO3A

nhood_markers5 <- findNhoodGroupMarkers(data_milo, data_res, subset.row = hvgs, 
                                        aggregate.samples = TRUE, sample_col = "orig.ident",
                                        subset.groups = c("5")
)

head(nhood_markers5)

#         logFC_5  adj.P.Val_5 GeneID
#MZB1    0.1608255 0.006357808   MZB1
#ANKUB1 -0.2244132 0.006357808 ANKUB1
#IGLC3   0.3006972 0.006357808  IGLC3
#NAA35  -0.3238196 0.006357808  NAA35
#DUSP1   0.3575025 0.006357808  DUSP1
#NFKBIA  0.2965150 0.006357808 NFKBIA

## We can do this to do differential expression between cells in different conditions within the same neighbourhood groups.
## only in those that have cells in both conditions (bad.prognosis | good prognosis) 

# - positive logFC means regulation
# - negative logFC means DOWNregulation

dge_2 <- testDiffExp(data_milo, data_res, design = ~ prognosis, meta.data = data.frame(colData(data_milo)),
                     subset.row = rownames(data_milo)[1:10], subset.nhoods=data_res$NhoodGroup=="2")
dge_2

dge_4 <- testDiffExp(data_milo, data_res, design = ~ prognosis, meta.data = data.frame(colData(data_milo)),
                     subset.row = rownames(data_milo)[1:10], subset.nhoods=data_res$NhoodGroup=="4")
dge_4

# we cannot test this in Nhood 5 as this neighbourhood has only cells from one condition (bad.prognosis)
# Add a column to distinguish the groups
nhood_markers5$group <- "Group 5"

# Rename columns to have consistent names across dataframes
colnames(nhood_markers5) <- c("logFC", "adj.P.Val", "GeneID", "group")

# Combine the dataframes

# Add a column to highlight significant genes based on adjusted p-value
nhood_markers5$Significant <- ifelse(nhood_markers5$adj.P.Val < 0.05, "Yes", "No")

# Identify the top 10 most significant genes
top10_genes <- nhood_markers5[order(nhood_markers5$adj.P.Val), ][1:50, ]

# Add a column for gene labels, with NA for non-top 10 genes
nhood_markers5$GeneLabel <- ifelse(nhood_markers5$GeneID %in% top10_genes$GeneID, nhood_markers5$GeneID, NA)

# Create a volcano plot
ggplot(nhood_markers5, aes(x = logFC, y = -log10(adj.P.Val), color = Significant, label = GeneLabel)) +
  geom_point() +
  geom_text_repel() +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Marker Genes", x = "Log Fold Change", y = "-log10 Adjusted P-Value") +
  theme(legend.title = element_blank()) +
  facet_wrap(~group)

# Filter the dataframe for significant genes (adj.P.Val < 0.05)
significant_genes <- combined_df[combined_df$adj.P.Val < 0.05, ]

# Extract and display the names of the significant genes
significant_gene_names <- significant_genes$GeneID
print(significant_gene_names)
