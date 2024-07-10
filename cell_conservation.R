
### Cell Cycle Conservation

cell_cycle_score <- function(seurat){
  
  # FUNCTION
  # The function computes the cell cycle score of each cell
  
  # INPUT
  # seurat <- seurat object
  
  # genes for evaluation of S- & G2M phases 
  s.genes <- cc.genes$s.genes # genes for S phase evaluation
  g2m.genes <- cc.genes$g2m.genes # genes for G2M evaluation 
  
  # scoring
  seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  
  return(seurat)
}


cell_cycle_conservation_score <- function(seurat_integrated, seurat_merged, mean = TRUE){
  
  # FUNCTION
  # The function computes a score evalutaing how well the cell-cycle effect can be captured 
  # before and after integration. 
  # It outputs a number between 0 (low conservation of variance explained by cell cycle)
  # and 1 (indicates )complete conservation of the variance).
  
  # INPUT
  # seurat_integrated <- integrated, processed seurat object
  # seurat_merged <- seurat object comprised by the merging and processing 
  #of the datasets used for integration
  # mean <- do you want the mean of the 2 values ?
  
  # score - 1-(var(cycle_score_integrated)-var(cycle_score_merged))/var(cycle_score_merged)
  
  # cell cycle score - S phase     
  var_s_integrated <- as.double(var(seurat_integrated[["S.Score"]]))
  var_s_merged <- as.double(var(seurat_merged[["S.Score"]]))
  cell_cycle_S_phase <- 1 - (abs(var_s_integrated - var_s_merged)/var_s_merged)
  
  # cell cycle score - G2M phase
  var_g2m_integrated <- as.double(var(seurat_integrated[["G2M.Score"]]))
  var_g2m_merged <- as.double(var(seurat_merged[["G2M.Score"]]))
  cell_cycle_G2M_phase <-  1 - (abs(var_g2m_integrated - var_g2m_merged)/var_g2m_merged)
  
  # plotting of relevant visualisation
  s_score <- FeaturePlot(seurat_integrated, features =  "S.Score") + labs(title = 'Feature plot of the S-phase scores in the integrated dataset')
  g2m_score <- FeaturePlot(seurat_integrated, features = "G2M.Score") + labs(title = 'Feature plot of the G2M-phase scores in the integrated dataset')
  print(s_score)
  print(g2m_score)
  
  if(mean == TRUE){
    
    # mean of scores
    cell_cycle_score <- (cell_cycle_S_phase + cell_cycle_G2M_phase)/2
    
    return(cell_cycle_score)
    
  }
  
  else{
    
    return(list(S_phase = cell_cycle_S_phase , G2M_phase = cell_cycle_G2M_phase))
  }
}
