make_chisq_heatmap <- function(x, reorder_values = NULL){
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(ComplexHeatmap)
  require(circlize)
  
  # Handle the matrices
  ph_res <- x$result[,c(1,2,4)] %>% pivot_wider(names_from = Cluster, values_from = res) %>% column_to_rownames("category") %>% as.matrix()
  ph_pva <- x$result[,c(1,2,3)] %>% pivot_wider(names_from = Cluster, values_from = pval) %>% column_to_rownames("category") %>% as.matrix()
  
  if (!is.null(reorder_values)){
    ph_res <- ph_res[,match(levels(reorder_values),colnames(ph_res))]
    ph_pva <- ph_pva[,match(levels(reorder_values),colnames(ph_pva))]
  }
  
  # stuff for heatmap
  col_enr_ramp <- 
    c(
      colorRampPalette(c(col_enr[1],col_enr[2]))(20)[1:19],
      colorRampPalette(c(col_enr[2],col_enr[3]))(21)
    )
  enr_pal_broad_1 <- circlize::colorRamp2(
    breaks = c(seq(min(ph_res),-0.1,length = 19),0,seq(0.1,max(ph_res),length = 20)), # this is NOT perfect. should be centered around zero no matter what. How to do this??
    colors = col_enr_ramp
  )
  pval_fun <- function(j,i,x,y,width,height,fill){
    if(ph_pva[i,j] < 0.01){
      if(ph_pva[i,j] < 0.001){
        grid.text("**", x, y, gp = gpar(fontsize = 10))
      } else{
        grid.text("*", x, y, gp = gpar(fontsize = 10))
      }
    }
  }
  
  # Make the heatmap
  ph_hm <-
    ComplexHeatmap::Heatmap(
      name = "residuals",
      ph_res,
      col = enr_pal_broad_1,
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      column_names_side = "top",
      column_names_gp = gpar(fontsize = 6),
      row_names_side = "left",
      cell_fun = pval_fun
    )
  
  # Gather output
  res <- list(
    residuals_matrix = ph_res,
    pvalues_matrix = ph_pva,
    heatmap = ph_hm
  )
  
  # return
  return(res)
}