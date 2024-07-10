### HVG (Highly Variable Genes) score

HVG_score <- function(dataset_integrated, dataset_merged, threshold = 500, assay_integrated = 'RNA', assay_merged = 'RNA'){
  
  
  # FUNCTION
  # The function computes the percentage of highy variable genes that are preserved 
  # through the integration process, here computed through the diffusion method. 
  # It is a proxy for the preservation of the biological signal.
  # Output is a number ranging from 0 (no highly variable genes are preserved through the process)
  # to 1 (all genes are preserved through the integration process).
  
  # INPUT
  # seurat_integrated <- integrated, processed seurat object
  # seurat_merged <- seurat object comprised by the merging and processing 
  #of the datasets used for integration
  # threshold <- how many genes do you want to keep for computation?
  # assay_integrated <- label of the assay in the integrated dataset for use
  # shouldn't coincide with assay_merged with integration techniques
  # that introduce another assay in the dataset
  # assay_merged <- label of the assay in the merged dataset
  
  # computation of features and gene expression
  # using "diffusion" method for higher applicability with different assays
  integrated_seurat_genes <- FindVariableFeatures(dataset_integrated, assay = assay_integrated, nfeatures = threshold, selection.method = 'dispersion')
  merged_seurat_genes <- FindVariableFeatures(dataset_merged, assay = assay_merged, nfeatures = threshold, selection.method = 'dispersion')
  
  # gene name extraction
  variable_test <- VariableFeatures(integrated_seurat_genes)
  variable_merged <- VariableFeatures(merged_seurat_genes)
  
  # scoring - len(intersect(genes_merged, genes_integrated))/min(len(genes_merged), genes_integrated)
  HVG_score <- length(intersect(variable_test, variable_merged))/(min(length(variable_merged), length(variable_test)))
  print(paste('The HVG conservation score is', as.character(HVG_score), sep = ' '))
  
  # Dot visualisation of 10 most present genes
  dot_integrated <- DotPlot(integrated_seurat_genes, features = VariableFeatures(integrated_seurat_genes)[1:10], group.by = "seurat_clusters") + ggtitle('DotPlot of the 10 most expressed genes in the integrated dataset, obtained through the diffusion method')
  dot_merged <- DotPlot(merged_seurat_genes, features = VariableFeatures(merged_seurat_genes)[1:10], group.by = "seurat_clusters") + ggtitle('DotPlot of the 10 most expressed genes in the merged dataset, obtained through the diffusion method')
  print(dot_integrated)
  print(dot_merged)
  
  return(HVG_score)
}