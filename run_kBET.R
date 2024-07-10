library(kBET)
library(Seurat)
library(ggplot2)

### kBET (k-nearest neighbour batch effect test)

run_kBET <- function(seurat, batch = 'Dataset', reduction = 'pca', downsample = TRUE, n_samples = 150){
  
  # FUNCTION
  # The function performs a batch effect test, evaluating the accordance of replicates 
  # based on Pearson's test chi-squared test for the label "batch" 
  # (which should only be a technical, non-biological grouping). 
  # The output is a model stored in a "list" object, with values rangeing from 0 to 1
  # and reppresent the likeliness that batch effects are present.
  
  # INPUT
  # seurat <- integrated, processed seurat object
  # batch <- batch for which to compute batch effect
  # should be a string of the name of the column in the metadata
  # downsample <- do you want to downsample ?
  # n_sample <- how many samples do you want to keep?
  
  
  # downsampling for home computer usability
  if (downsample){
    seurat <- subset(seurat, downsample = n_samples) 
  }
  # case use - only 1 element in the selected label
  if(length(unique(unlist(seurat[[batch]])))<2){
    print(paste('The dataset contains only label', unique(unlist(seurat[[batch]])), 'for batch', batch, '- computation not performed', sep = ' '))
    return(NA)
  }
  
  # extraction of reduction embeddings and of batch variable
  print('na')
  pca_embeddings <- Embeddings(seurat, reduction = "pca")
  print('ba')
  batch_variable <- unlist(seurat[[batch]])
  
  # kBET computation
  print('ja')
  kBET_results <- kBET(pca_embeddings, batch_variable, plot = FALSE)
  print('bo')
  print(kBET_results$summary)
  print('ggg')
  
  # extraction & plotting of results
  plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                    each=length(kBET_results$stats$kBET.observed)), 
                          data =  c(kBET_results$stats$kBET.observed,
                                    kBET_results$stats$kBET.expected))
  kbet_plot <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
    labs(x='Test', y='Rejection rate',title=paste('kBET test results for batch',
                                                  batch, sep = ' ')) + theme_bw() + scale_y_continuous(limits=c(0,1))
  print(kbet_plot)
  
  return(kBET_results$summary)
}
