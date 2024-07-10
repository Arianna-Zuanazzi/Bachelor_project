library(cluster)

### ASW ( Average Width Silhouette) & F1

get_isolated_labels <- function(seurat, celltype_label = "celltype_major", sample_label = "orig.ident", iso_threshold = 5) {
  
  
  # FUNCTION
  # The function extracts and outputes the isolated labels from a dataset, 
  # defined as the labels present in less than iso_threshold samples, 
  # in the form of a vector of strings.
  
  # INPUT
  # seurat <- integrated, processed seurat object
  # celltype_label <- label of the celltype
  # should be a string with the name of the column in the metadata
  # sample_label <- label of the sample
  # should be a string with the name of the column in the metadata
  # iso_threshold <- in how many samples max should the celltype be ?
  
  # extraction of celltype labels per sample
  temp <- unique(seurat@meta.data[, c(celltype_label, sample_label)])
  batch_per_lab <- aggregate(temp[, sample_label], by = list(temp[, celltype_label]), FUN = length)
  colnames(batch_per_lab) <- c(celltype_label, sample_label)
  #threshold for the label to be considered isolated
  if (is.null(iso_threshold)) {
    iso_threshold <- min(batch_per_lab[, sample_label])
  }
  
  # extraction of true isolated labels 
  labels <- batch_per_lab[batch_per_lab[, sample_label] <= iso_threshold, celltype_label]
  if (length(labels) == 0) {
    print(paste0("no isolated labels with less than ", iso_threshold, " batches\n"))
  }
  
  return(labels)
}


isolated_labels_dataset <- function(seurat, labels){
  
  
  # FUNCTION
  # The function subsets an integrated seurat object
  # leaving in only the celltypes present in the isolated labels.
  
  # INPUT
  # seurat <- integrated, processed seurat object
  # labels <- celltypes present in a scarce, predefined number of samples
  # should be a vector of strings with the names of the celltypes to keep
  
  # creation of isolated labels dataset - keeping only cells with celltype in labels
  isolated_dataset <- subset(seurat, celltype_major %in% labels)
  
  return(isolated_dataset)
}

ASW_score <- function(seurat, batch = 'Dataset', reduction = 'pca', downsample = TRUE, n_samples = 800){
  
  
  # FUNCTION
  # The function computed the average silhouette width of a dataset for a specific label.
  # Using biological labels (celltype, cancer type, ...)
  # it gives a measure of the conservation of biological variance, 
  # using technical labels (dataset, sequencing technology ...)
  # it gives a measure of the removal of batch effects. 
  # It returns a number ranging from -1 (high dissimilarity between elements of the same cluster)
  # up to 1 (high similarity).
  
  # INPUT
  # seurat <- integrated, processed seurat object
  # batch <- batch for which to compute ASW
  # should be a string with the name of the column in the metadata
  # reduction <- name of the reduction chosen
  # should be a string with the name of the reduction (usually pca, umap or tsne)
  # downsample <- do you want to downsample ?
  # n_samples <- down to how many samples ?
  
  # downsampling for home computer usability
  if(downsample){
    seurat <- subset(seurat, downsample = n_samples)
  }
  
  # cause use - only 1 label in selected batch 
  if(length(unique(unlist(seurat[[batch]])))== 1){
    print(paste('ASW not computable for label', batch,'for only one unique instance of the label is present in the dataset', sep = ' '))
    return(NA)
  }
  
  # computation of dissimilarity matrix from reduction embeddings
  dist.matrix <- dist(x = Embeddings(object = seurat[[reduction]]))
  # extraction of batch label & silhouette computation
  batch_labels <- as.numeric(as.factor(unlist(seurat[[batch]])))
  silhouette <- silhouette(x = batch_labels, dist = dist.matrix)
  
  # ASW computation - mean of silhouette values
  ASW <- silhouette[, 3] %>% mean(na.rm = TRUE)
  print(paste("The ASW score for the categorical label", batch, "is", as.character(ASW), sep = ' '))
  
  # plotting of silhouette
  sil <- fviz_silhouette(silhouette)+ scale_fill_discrete(labels = unique(unlist(seurat[[batch]]))) + ggtitle(paste('Silhouette of the dataset, clustered by', batch, sep = ' '))
  print(sil)
  
  return(ASW)
}


labels_F1 <- function(seurat, batch_prec = 'celltype_major', batch_rec = 'seurat_clusters', labels = c()){
  
  # FUNCTION
  # The function computes measure of a test's accuracy, 
  # computed as the harmonic mean between precision and recall. 
  # It returns a value ranging from 0 to 1, 
  # indicating how often the test is making the right prediction.
  
  
  # INPUT
  # seurat <- integrated, processed seurat object
  # batch_prec <- batch used as the "precision" factor, usually celltype
  # should be a string with the name of the column in the metadata
  # batch_prec <- batch used as the "recall"/"sensitivity" factor, usually cluster
  # should be a string with the name of the column in the metadata

  # computaton of normalised precision & recall ASW 
  celltype_ASW <- (ASW_score(seurat, batch = batch_prec)+1)/2
  cluster_ASW <- (ASW_score(seurat, batch = batch_rec)+1)/2
  
  # computation of F1 score
  F1 <- 2*celltype_ASW*cluster_ASW/(celltype_ASW + cluster_ASW)
  print(paste('The isolated label F1 is', as.character(F1), sep = ' '))
  
  
  return(F1)
}
