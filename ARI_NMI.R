library(aricode)

### ARI and NMI (Adjusted Rand Index & Normalised Mutual Information)

ARI_batch <- function(seurat, batch_label = 'Dataset', cluster_label = 'seurat_clusters'){
  
  # FUNCTION
  # The function indicatively computes the measure of the similarity between two data clusterings,
  # adjusted by the possibility that it may arise by chance. 
  # It outputes a number ranging from 0 (-0.5 for especially bad correlation) to 1
  # showing the degree of degree of overlapping between labels
  
  # INPUT
  # seurat <- integrated, processed seurat object
  # batch_label <- batch for ARI computation
  # should be a string with the name of the column in the metadata
  # cluster_label <- clustering of reference
  # should be a string with the name of the column in the metadata
  
  
  # ARI computation
  ari_result <- ARI(unlist(seurat[[batch_label]]), unlist(seurat[[cluster_label]]))
  print(paste('The ARI score between batch', batch_label, 'and batch', cluster_label, 'is', as.character(ari_result), sep = ' '))
  
  return(ari_result)
}

NMI_batch <- function(seurat, batch_label = 'Dataset', cluster_label = 'seurat_clusters'){
  
  
  # FUNCTION
  # The function indicatively computes the quality of clustering compared to another, 
  # and accounts to the “amount of information” one can extract 
  # from a distribution regarding a second one.
  # It outputes a number ranging from 0 (no mutual information) and 1 (perfect correlation).
  
  # INPUT
  # seurat <- integrated, processed seurat object
  # batch_label <- batch for ARI computation
  # should be a string with the name of the column in the metadata
  # cluster_label <- clustering of reference
  # should be a string with the name of the column in the metadata
  
  
  # NMI computation
  nmi_result <- NMI(unlist(seurat[[batch_label]]), unlist(seurat[[cluster_label]]))
  print(paste('The NMI score between batch', batch_label, 'and batch', cluster_label, 'is', as.character(nmi_result), sep = ' '))
  
  return(nmi_result)
}
