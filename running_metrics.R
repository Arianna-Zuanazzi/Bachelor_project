# RUNNING METRICS

integrated_datasets <- c(rpca_integrated, fastmnn_integrated, harmony_integrated, liger_integrated)
datasets <- append(integrated_datasets, processed_merged_dataset)

# Processing for metrics
for(dataset in datasets){
  dataset <- pseudotime_monocle3(dataset)
  dataset <- pseudotime_pca(dataset)
  dataset <- cell_cycle_score(dataset)
}

for(dataset in integrated_datasets){
  dataset_name <- sub('.*_', '', levels(Idents(dataset)))
  print(paste('Now processing dataset', dataset_name, sep = ' '))
  print(' ')
  name <- paste('Dataset', dataset_name, 'metrics', sep = '_')
  metrics <- metrics_panel(dataset, processed_merged_dataset)
  assign(name, metrics)
}
