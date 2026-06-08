#quick featureplot
feature <- function(x,scdata){
  FeaturePlot(
    scdata,
    reduction = "umap",
    dims = c(1,2),
    features = x,
    cols = c("lightgrey","#a100d2"),
    pt.size = 0.75,
    order = T,
    label = TRUE,
    slot = "counts",
    label.size = 2.5
  )
}
