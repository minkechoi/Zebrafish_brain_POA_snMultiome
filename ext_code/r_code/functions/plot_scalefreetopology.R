plot_scalefreetopology <- function(sft){
  par(mfrow = c(1, 2))
  cex1 = 0.9
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    ylim = c(0,1),
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit,signed R^2",
    pch = ".",
    main = paste("Scale independence")
  )
  
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    pos = 3,
    labels = powers,
    cex = cex1,
    col = "darkred"
  )
  
  # this line corresponds to using an R^2 cut-off of h
  abline(h = 0.85, col = "tomato")
  abline(h = 0.9, col = "lightgreen")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    pch = ".",
    main = paste("Mean connectivity")
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    pos = 3,
    cex = cex1,
    col = "darkred"
  )
  par(mfrow = c(1,1))
}


plot_scalefreetopology_pretty <- function(sft){
  require(colorspace)
  require(circlize)
  give_col = circlize::colorRamp2(
      breaks = quantile(c(-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2])),
      col = sequential_hcl(5, "YlGnBu", rev = TRUE)
    )

  par(mfrow = c(1, 2))
  cex1 = 0.9
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    ylim = c(0,1),
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit,signed R^2",
    pch = 21,
    type = "b",
    bg = give_col(-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]),
    bty = "n",
    main = paste("Scale independence")
  )
  
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    pos = 3,
    labels = powers,
    cex = cex1
  )
  
  # this line corresponds to using an R^2 cut-off of h
  abline(h = 0.85, col = "tomato", lwd = 1.2)
  abline(h = 0.9, col = "mediumseagreen", lwd = 1.2)
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    pch = 21,
    type = "b",
    bg = "gray95",
    col = darken("gray95",.5),
    bty = "n",
    main = paste("Mean connectivity")
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    pos = 3,
    cex = cex1
  )
  par(mfrow = c(1,1))
}
