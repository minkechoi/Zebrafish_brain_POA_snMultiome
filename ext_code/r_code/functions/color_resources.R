#' Color and graphic materials (functions, palettes, resources, ...)

# Color palettes
viridis_pastel <- c("#ffee61", "#96e88e", "#5dc9ac","#4da2ba","#6b6eab","#552761")

# from https://stackoverflow.com/questions/50600425/r-convert-colors-to-pastel-colors
pastelise_hsv <- function(col,n=0.4){
  y = hsv(col[1], col[2]*n, col[3])
  return(y)
}
# Works in loop functions context. Example of how to use:
# vector_palette = c("red","blue","green","yellow")
# sapply(as.data.frame(rgb2hsv(col2rgb(vector_palette))),pastelise_hsv)

gb_col <- 
  c(
    yellow1 = "#f5f8aa",
    yellow2 = "#d0be7b",
    orange1 = "#efad84",
    orange2 = "#d78c91",
    purple1 = "#9d95e2",
    purple2 = "#7b6b94",
    blue1 = "#92baef",
    blue2 = "#4a42b5",
    teal1 = "#6bad94",
    teal2 = "#2f738a",
    green1 = "#9cce7b",
    green2 = "#4a9c5a",
    white1 = "#eff0f9",
    gray1 = "#9696a4",
    gray2 = "#5a5a63",
    black1 = "#1a1a1a"
  )

average_cols <- function(x){ # from https://stackoverflow.com/a/29576746
  if(length(x)<2) stop("only one color.")
  a <- col2rgb(x)
  b <- sqrt(sapply(data.frame(t(a)),function(x){mean(x^2)}))/255
  y <- rgb(b[1],b[2],b[3])
  return(y)
}
