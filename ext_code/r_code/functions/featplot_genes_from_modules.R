
#' @params module_list: string of list of modules
#' @params scdata: single cell object to plot the genes
#' @params output_rooth_path: the (pre-existing) path where the rest of the folders will be created
#' @params num_genes: minimum number of genes to plot. Defaults to 30.
#' @params size: png size. Default is a square same-sized png file. Defaults to 256 pixels
#' @params random_seed: the random seed to use for subsampling the genes to plot per module
featplot_genes_from_modules <- 
  function(
    module_list, wg_module, scdata, output_root_path,
    num_genes=30, size = 256, random_seed = 5678
    ){
    
    require(Seurat)
    
    for(i in module_list){
      
      # Create directory
      newfolder <- paste0("module_",i)
      newpath <- file.path(output_root_path, newfolder)
      dir.create(newpath)
      
      
      set.seed(random_seed)
      all_i <- rownames(wg_module[wg_module$module == i,])
      
      if(length(all_i) < num_genes ){
        g_i <- all_i
      } else {
        g_i <- sample(all_i, num_genes)
      }
      
      print(i)
      
      for(j in g_i){
        
        j_n <- which(g_i == j)
        
        filename <-
          paste0(
            newpath,"/",
            "umap_module_",
            i,
            "_",
            j_n,
            ".png"
          )
        
        featplot_g_i <- 
          FeaturePlot(
            scdata,
            features = j,
            order = TRUE,
            raster = FALSE,
            pt.size = 1
          ) + NoLegend() + NoAxes()
        
        png(
          file = filename,
          width= size,
          height = size
        )
        print(featplot_g_i)
        dev.off()
        
        message("done ",j)
        
      }
      
      message("done ", i)
      
    }
  }
