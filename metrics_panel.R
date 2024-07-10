metrics_panel <- function(seurat_integrated, seurat_merged){
  
  # FUNCTION
  # The function performs a full panel of metrics testing.
  # It asks for the names of the following label inputs (as strings)
  # dataset AUTHOR
  # SAMPLE within the dataset
  # CELLTYPE
  # CLUSTER
  # single cell RNA SEQUENCING PROCESS
  # SEQUENCING TECHNOLOGY 
  # DIMENSIONALITY REDUCTION (usually pca, umap or tsne)
  # ASSAY of the INTEGRATED dataset 
  # ASSAY of the MERGED dataset (usually RNA)
  # PSEUDOTIME
  # whether STRICT and LASS DOWNSAMPLING should be performed
  # and down to how many cells
  # PRECISION and SENSITIVITY labels
  # Different labels will be used on different metrics, divided in groups
  # Technical labels = author, sample, sequencing technology, sequencing RNA process
  # Batches = technical labels + celltype
  # Labels = Batches + cluster
  # The following metrics are implemented
  # Principal Component Regression (on batches)
  # kBET (on technical batches)
  # ARI (on batches, comparing label = cluster)
  # NMI (on batches, comparing label = cluster)
  # ASW (on labels)
  # isolated label ASW (on labels)
  # LISI (on labels)
  # F1 (on celltype + cluster)
  # Trajectory conservation
  # HVG score
  # It requires for the two datasets to be seurat objects, filtered and processed as for standard
  # with no data scaling performed. Pseudotime and cell cycle scoring need to have been performed, 
  # as well as PCA and UMAP reduction.
  # Output is a nested list with all metrics uniquely named and classified according to macro type 
  
  # INPUT
  # seurat_integrated <- integrated, processed seurat object
  # seurat_merged <- seurat object comprised by the merging and processing 
  #of the datasets used for integration
  
  
  # establishment of metrics
  print('Establishment of required variables for metric purposes')
  
  print('Write the name of the AUTHOR label.')
  author <- readline('If you leave an empty line, the standard labels "authors" will be used: ')
  if(author == ''){author <- 'authors'}
  
  
  print('Write the name of the SAMPLE label.')
  sample <- readline('If you leave an empty line, the standard labels "orig.ident" will be used: ')
  if(sample == ''){sample <- 'orig.ident'}
  
  print('Write the name of the CELLTYPE label.')
  celltype <- readline('If you leave an empty line, the standard label "celltype_major" will be used: ')
  if(celltype == ''){celltype <- 'celltype_major'}
  
  
  print('Write the name of the CLUSTER label.')
  cluster <- readline('If you leave an empty line, the standard labels "seurat_clusters" will be used: ')
  if(cluster == ''){cluster <- 'seurat_clusters'}
  
  print('Write the name of the SINGLE CELL RNA TECHNOLOGY label.')
  seq_type <- readline('If you leave an empty line, the standard labels "seq_type" will be used: ')
  if(seq_type == ''){seq_type <- 'seq_type'}
  
  print('Write the name of the PROCESSING TECHNOLOGY label.')
  seq_tech <- readline('If you leave an empty line, the standard labels "seq_tech" will be used: ')
  if(seq_tech == ''){seq_tech <- 'seq_tech'}
  
  print('Write the name of the REDUCTION label.')
  reduction <- readline('If you leave an empty line, the standard labels "pca" will be used: ')
  if(reduction == ''){reduction <- 'pca'}
  
  print('Write the name of the ASSAY of the INTEGRATED dataset.')
  assay_integrated <- readline('If you leave an empty line, the standard label "RNA" will be used: ')
  if(assay_integrated == ''){assay_integrated <- 'RNA'}
  
  
  print('Write the name of the ASSAY of the MERGED dataset.')
  assay_merged <- readline('If you leave an empty line, the standard label "RNA" will be used: ')
  if(assay_merged == ''){assay_merged <- 'RNA'}
  
  
  print('Write the name of the PSEUDOTIME label.')
  pseudotime <- readline('If you leave an empty line, the standard label "pseudotime_monocle" will be used: ')
  if(pseudotime == ''){pseudotime <- 'pseudotime_monocle'}
  
  
  
  large_downsample <- readline('Write "y" if you do want to perform LARGE DOWNSAMPLING for the metrics ASW and LISI:')
  if(large_downsample %in% c('Y', 'y')){large_downsample <- TRUE}
  
  
  if(large_downsample == TRUE){
    print('Write down numerically how many samples do you want to downsample to.')
    large_samples <- readline('If you leave an empty line, the standard number 800 will be used: ')
    if(large_samples == ''){large_samples <- 800}
    else{large_samples <- as.numeric(large_samples)}
  }
  
  
  strict_downsample <- readline('Write "y" if you do want to perform STRICT DOWNSAMPLING for the metrics kBET:')
  if(strict_downsample %in% c('y', 'Y')){strict_downsample <- TRUE}
  
  
  if(strict_downsample == TRUE){
    print('Write down numerically how many samples do you want to downsample to.')
    strict_samples <- readline('If you leave an empty line, the standard number 150 will be used: ')
    if(strict_samples == ''){strict_samples <- 150}
    else{strict_samples <- as.numeric(strict_samples)}
  }
  
  print('Write down numerically the GENE THRESHOLD you wan to use for HVG score computation.')
  gene_threshold <- readline('If you leave an empty line, the standard number 500 will be used: ')
  if(gene_threshold == ''){gene_threshold <- 500}
  else{gene_threshold <- as.numeric(gene_threshold)}
  
  print('Write down numerically the THRESHOLD you want to use to determine ISOLATED LABELS.')
  isolated_labels_threshold <- readline('If you leave an empty line, the standard number 5 will be used: ')
  if(isolated_labels_threshold == ''){isolated_labels_threshold <- 5}
  else{isolated_labels_threshold <- as.numeric(isolated_labels_threshold)}
  
  print('Write down the name of the PRECISION label to be used for the computing of the isolated lables F1 score.')
  precision_label <- readline('If you leave an empty line, the celltype label will be used: ')
  if(precision_label == ''){precision_label <- celltype}
  
  print('Write down the name of the RECALL label to be used for the computing of the isolated lables F1 score.')
  recall_label <- readline('If you leave an empty line, the cluster label will be used: ')
  if(recall_label == ''){recall_label <- cluster}
  
  
  
  # checking & printing variables
  variables <- c(author, sample, celltype, cluster, reduction, assay_integrated, assay_merged, pseudotime, as.character(large_downsample), as.character(strict_downsample), as.character(gene_threshold), as.character(isolated_labels_threshold), precision_label, recall_label, seq_type, seq_tech)
  print('These are the names of the variables that are gonna be used')
  lapply(variables, print)
  
  # group establishment - technical batches, general labels for comparison, all labels
  technical_batches <- c(author, sample, seq_tech, seq_type)
  batches <- c(author, sample, seq_tech, seq_type, celltype)
  labels <- c(author, sample, celltype, seq_tech, seq_type, cluster)
  
  # metrics variable
  metrics <- list()
  
  # PCR
  pcr_model <- run_pcr(seurat_integrated, batches = batches, target = cluster)
  metrics[['PCR_model']] <- pcr_model
  
  # kBET
  kBET_scores <- list()
  for(batch in c(author, seq_tech, seq_type)){
    name <- paste('kBET', batch, sep = '_')
    score <- run_kBET(seurat_integrated, batch = batch, reduction = reduction, 
                      downsample = strict_downsample, n_samples = strict_samples)
    kBET_scores[[name]] <- score
  }
  metrics[['kBET_scores']] <- kBET_scores
  
  
  # ARI
  ARI_scores <- list()
  for(batch in batches){
    x <- ARI_batch(seurat_integrated, batch_label = batch, cluster_label = cluster)
    name <- paste('ARI', batch, sep = '_')
    ARI_scores[[name]] <- x
  }
  metrics[['ARI_scores']] <- ARI_scores
  
  # NMI
  NMI_scores <- list()
  for(batch in batches){
    score <- NMI_batch(seurat_integrated, batch_label = batch, cluster_label = cluster)
    name <- paste('NMI', batch, sep = '_')
    NMI_scores[[name]] <- score
  }
  metrics[['NMI_scores']] <- NMI_scores
  
  
  # LISI
  LISI_scores <- list()
  for(label in labels){
    score <- LISI_score(seurat_integrated, reduction_type = reduction, 
                        label = label, downsample = large_downsample, n_samples = large_samples)
    name <- paste('LISI', label, sep = '_')
    LISI_scores[[name]]<- score
  }
  metrics[['LISI_scores']] <- LISI_scores
  
  
  # ASW
  ASW_scores <- list()
  
  # ASW, non-isolated labels
  ASW_non_isolated <- list()
  for(label in labels){
    name <- paste('ASW', label, spe = '_')
    score <- ASW_score(seurat_integrated, batch = label, reduction = reduction, downsample = large_downsample, n_samples = large_samples)
    ASW_non_isolated[[name]] <- score
  }
  ASW_scores[['ASW_non_isolated']] <- ASW_non_isolated
  
  # ASW, isolated
  ASW_isolated <- list()
  isolated_labels <- get_isolated_labels(seurat_integrated, celltype_label = celltype, sample_label = sample, 
                                         iso_threshold = isolated_labels_threshold)
  isolated_dataset <- isolated_labels_dataset(seurat_integrated, labels = isolated_labels)
  for(label in labels){
    name <- paste('ASW', 'isolated',label, sep = '_')
    score <- ASW_score(isolated_dataset, batch = label, reduction = reduction, downsample = large_downsample, n_samples = large_samples)
    ASW_isolated[[name]] <- score
  }
  
  # ASW, isolated F1
  name <- 'ASW_F1_isolated'
  score <- isolated_labels_F1(isolated_dataset, batch_prec = precision_label, batch_rec = recall_label)
  ASW_isolated[name] <- score
  ASW_scores[['ASW_isolated']] <- ASW_isolated
  metrics[['ASW_scores']] <- ASW_scores
  
  
  # Trajectory conservation
  name <- 'trajectory_conservation'
  score <- trajectory_conservation(seurat_integrated, seurat_merged, pseudotime_label = pseudotime)
  metrics[[name]] <- score
  
  # Cell Cycle conservation
  name <- 'cell_cycle_conservation'
  score <- cell_cycle_conservation_score(seurat_integrated, seurat_merged, mean = TRUE)
  metrics[['cell_cycle_score']] <- score
  
  # HVG score
  name <- 'HVG_score'
  score <- HVG_score(seurat_integrated, seurat_merged, threshold = gene_threshold, 
                     assay_integrated = assay_integrated, assay_merged = assay_merged)
  metrics[[name]] <- score
  
  return(metrics)
}