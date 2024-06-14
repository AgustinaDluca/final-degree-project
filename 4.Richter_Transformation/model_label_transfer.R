
library(matchSCore2)
library(nnet)
library(Matrix)

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

## To each cell probabilities are assigned for any possible identity class
probabilities <- out$fit.prob



