library(pls)

### Principal Component Regression

run_pcr <- function(seurat, batches = c('authors', 'orig.ident', 'celltype_major'), target = 'seurat_clusters'){
  
  # FUNCTION
  # The function performs a principal component regression of the dataset out of the batches specified
  # under "batches", predicting vales in the variables stored under "target",
  # returning a model stored as a "list" object.
  # The results of this function can be used to assess number of PCA to keep
  # as well as predict and fit future data (using function predict() ).
  
  # INPUT
  # seurat <- integrated& processed seurat object
  # batches <- batches for which to compute pcr
  # should be a vector containing the names of the columns in the metadata for the pcr computation
  # target <- column in the metadata for which to compute pcr
  # should be a string of the column name in the metadata
  
  
  # dataframe extraction, subsetting to relevant columns
  dataframe <- as.data.frame(as.matrix(seurat@meta.data))
  columns <- append(batches, target)
  dataframe <- dataframe[, columns]
  
  # elements of column turned into numeric values, extraction of target variable
  for(batch in batches){dataframe[[batch]] <- as.numeric(as.factor(dataframe[[batch]]))}
  target_variable <- dataframe[[target]]
  
  # building of pcr model, plotting relevant visualisations
  pcr_model <- pcr(target_variable ~ ., data = dataframe, validation = "CV",)
  print(summary(pcr_model))
  print(validationplot(pcr_model, val.type = "MSEP"))
  
  return(pcr_model)
}
