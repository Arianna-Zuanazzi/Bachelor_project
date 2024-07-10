### Pseudotime & Trajectory Conservation

pseudotime_monocle3 <- function(dataset, num_dimension = 40){
  
  
  # FUNCTION
  # The functions uses Monocle3 to computed pseudotime, 
  # the latent (unobserved) dimension which measures the cells’ progress through the transition.
  # It returns a seurat object containing the pseudotime value (numerical values) for each cell.
  
  # INPUT
  # dataset <- dataset to be used, prefereably a seurat object
  # num_dimensions <- number of dimensions used in the internal monocle3 process 
  
  # consersion of seurat object into cell dataset 
  cds <- as.cell_data_set(dataset)
  
  # preprocessiong, clustering and graphing
  cds <- preprocess_cds(cds, num_dim = num_dimension)
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  
  # computing pseudotime values
  cds <- order_cells(cds)
  # Adding pseudotime values back to Seurat object
  dataset[["monocle_pseudotime"]] <- cds@principal_graph_aux[["UMAP"]]$pseudotime
  
  return(dataset)
}

pseudotime_pca <- function(seurat){
  
  
  # FUNCTION
  # The functions assigns pseudotimes, the latent (unobserved) dimension 
  # which measures the cells’ progress through the transition,
  # definined as the values of the 1st Principal Component for the cell.
  # It returns a seurat object containing the pseudotime value (numerical values) for each cell.
  
  # INPUT
  # seurat <- integrated, processed seurat object
  
  # extraction of PC_! values from cells & adding to seurat object 
  pca_pseudotime <- FetchData(seurat, 'PC_1')
  seurat[['pseudotime_pca']] <- pca_pseudotime
  return(seurat)
}

trajectory_conservation <- function(seurat_integrated, seurat_merged, pseudotime_label = 'monocle_pseudotime'){
  
  # FUNCTION
  # The function computes a proxy value for the conservation of the biological signal. 
  # It returns a number between 1 (the cells lie on the same order on the trajectory 
  # before and after integration)
  # and 0 (the cells lie on reverse order on the trajectory).
  
  # INPUT
  # seurat_integrated <- integrated, processed seurat object
  # seurat_merged <- seurat object comprised by the merging and processing 
  #of the datasets used for integration
  # pseudotime_label <- label of the pseudotime values
  # should be a string with the name of the column in the metadata
  
  # extracting pseudotime values 
  pseudotime_integrated <- unlist(seurat_integrated[[pseudotime_label]])
  pseudotime_merged <- unlist(seurat_merged[[pseudotime_label]])
  
  # computing correlation & extraction of the correlation coefficient
  result <- cor.test(pseudotime_integrated, pseudotime_merged, method = "spearman")
  rho <- result$estimate
  
  # scaling the score to a value between 0 and 1
  trajectory_conservation_score <- (rho + 1) / 2
  print(paste('The trajectory conservation score is', as.character(trajectory_conservation_score), sep = ' '))
  
  # plotting pseudotime values
  pseudo_integrated <- FeaturePlot(seurat_integrated, features = pseudotime_label) + labs(title = 'FeaturePlot of the pseudotime in the integrated dataset')
  pseudo_merged <- FeaturePlot(seurat_merged, features = pseudotime_label) + labs(title = 'FeaturePlot of the pseudotime in the merged dataset')
  print(pseudo_integrated)
  print(pseudo_merged)
  
  return(trajectory_conservation_score)
}