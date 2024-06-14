# BiocManager::install("UCell")
library(UCell)

signature <- list("CLL_bad.prognosis_signature"= c("ZAP70", "CD38","CD27","CXCR4","IGHM", "MZB1","XBP1", 'HDAC9', 'HSD17B12', 'GTF2E2'))

data <- AddModuleScore_UCell(data, features = signature, name=NULL)
FeaturePlot(data_rt, reduction = "umap.rna", features = names(signature), min.cutoff = 0.15, cols = cols_sign, order = F)

VlnPlot(data, group.by = 'new_annotation', features = names(signature), sort = T)

# Create the dot plot to get the average expression values
plot <- DotPlot(data1, features = names(signature), group.by = 'new_annotation') + RotatedAxis() + coord_flip()
plot_data <- plot$data
avg_exp <- tapply(plot_data$avg.exp.scaled, plot_data$id, mean)
ordered_clusters <- factor(names(avg_exp), levels = names(sort(avg_exp)))

data1@meta.data$new_annotation <- factor(data1@meta.data$new_annotation, levels = levels(ordered_clusters))

a = DotPlot(data, features = names(signature), group.by = 'new_annotation', cols = cols_sign, dot.scale = 10, col.min = -1.5, col.max = 0.5) + RotatedAxis() + coord_flip()
a = as.data.frame(a$data)

ggplot(a, aes(x = features.plot, y = id, fill = avg.exp)) +
  geom_tile(color = "gray") +  # Add white borders around tiles for better visibility
  scale_fill_gradient(low = "yellow1", high = "purple", limits = c(0.18, 0.32), na.value = "grey50", name = "Average Expression") +  # Adjust color gradient limits
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + coord_flip() + RotatedAxis()

b = DotPlot(data, features = signature, group.by = 'new_annotation', cols = cols_sign, dot.scale = 10, col.min = -1.5, col.max = 0.5) + RotatedAxis() + coord_flip()
b = as.data.frame(b$data)

ggplot(b, aes(x = features.plot, y = id, fill = avg.exp.scaled)) +
  geom_tile(color = "gray") +  # Add white borders around tiles for better visibility
  scale_fill_gradient(low = "yellow", high = "purple", limits = c(-1.5, 0.5), na.value = 'gray', name = "Average Expression") +  # Adjust color gradient limits
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + RotatedAxis()
