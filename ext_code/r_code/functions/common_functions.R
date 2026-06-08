load_matrices <- function(tablemtx, min_cells = 1, min_feats = 100, ...){
  
  library(vroom)
  library(Seurat)
  library(tibble)
  
  seur = list()
  
  colnames(tablemtx)[1:4] <- c("unit","project","sublibrary","path")
  
  tablemtx$mtx_units <- paste0(
    tablemtx$unit,
    "_",
    tablemtx$project,
    ".",
    tablemtx$sublibrary
  )
  
  for (i in 1:nrow(tablemtx)) {
    print(paste0("Loading ",tablemtx$mtx_units[i],"..."))
    
    mtx <-
      vroom(
        file = tablemtx$path[i]
      )
    
    mtx <- mtx %>% column_to_rownames(var = "GENE")
    
    print(paste0("Creating SeuratObject ",tablemtx$mtx_units[i],"..."))
    
    seur[i] <- 
      CreateSeuratObject(
        counts = mtx,
        project = tablemtx$mtx_units[i],
        min.cells = min_cells,
        min.features = min_feats
        )
      
    names(seur)[i] <- tablemtx$mtx_units[i]
    
    }
  
  return(seur)
  
}


fcha <- function(){ gsub("-","",Sys.Date()) }