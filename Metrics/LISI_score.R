library(lisi)


# LISI (Local Inverse Simpson's Index)

LISI_score <- function(seurat, reduction_type = 'pca', label = 'authors', downsample = TRUE, n_samples = 800){
  
  
  # FUNCTION
  # The function assess degree of mixing during batch correction and dataset integration
  # and computes the effective number of different categories 
  # represented in the local neighborhood of each cell.
  # Using biological labels (celltype, cancer type, ...)
  # it gives a measure of the conservation of biological variance, 
  # using technical labels (dataset, sequencing technology ...)
  # it gives a measure of the removal of batch effects. 
  # It outputs a number ranging to 1 to the number of possible elemnts present in the chosen batch.
  
  # INPUT
  # seurat <- integrated, processed seurat object
  # reduction_type <- reduction chosen
  # should be a string with the exact name of the reduction chosen (usually pca, umap or tsne)
  # label <- label for which to compute LISI
  # should be a string with the name of the column in the metadata
  # downsample <- do you wanna downsample ?
  # n_samples <- down to how many samples ?
  
  # downsampling for home computer usability   
  if(downsample){
    seurat <- subset(seurat, downsample = n_samples)
  }
  
  # extracion of embeddings
  reduction_embeddings <- Embeddings(object = seurat, reduction = reduction_type)
  
  # computation of lisi values 
  lisi_values <- compute_lisi(reduction_embeddings, seurat@meta.data, label_colnames = label)
  
  # LISI computation, mean of lisi values 
  lisi_metric <- mean(lisi_values[, 1])
  print(paste('The LISI score for the categorical label', as.character(label), 'is', as.character(lisi_metric), sep = ' '))
  
  return(lisi_metric)
}
