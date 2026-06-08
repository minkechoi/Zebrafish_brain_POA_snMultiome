# Create Matrix
set.seed(1234)
mat <- cbind(
  data.frame(
    ATAC1 = rpois(150, lambda = 5),
    ATAC2 = 2+rpois(150, lambda = 5),
    ATAC3 = 3+rpois(150, lambda = 5),
    ATAC4 = rpois(150, lambda = 5),
    ATAC5 = rpois(150, lambda = 5)
  ),
  data.frame(
    RNA1 = rpois(150, lambda = 5),
    RNA2 = 3+rpois(150, lambda = 5),
    RNA3 = 4+rpois(150, lambda = 5),
    RNA4 = rpois(150, lambda = 5),
    RNA5 = rpois(150, lambda = 5)
  )
)
# Create Matrix; add negative_correlation
set.seed(1234)
sss = sample(1:nrow(mat), 50)
mat[sss,c(2,3)] <- abs(mat[sss,c(2,3)]-5)
mat[sss,c(1,4)] <- abs(mat[sss,c(2,3)]+5)
  
# heatmap
h1 <- Heatmap(name = "FC", t(scale(t(mat))), cluster_columns = FALSE, show_row_names = FALSE, column_split = c("ATAC","ATAC","ATAC","ATAC","ATAC","RNA","RNA","RNA","RNA","RNA"), column_title = "Pos & Neg correlation")

# correlations between atac and rna
cor_atac_rna <-
  apply(
    mat,
    1,
    function(x){
      cor(x[1:5],x[6:10])
    }
  )
# heatmap only with positive correlation
h2 <- Heatmap(name = "FC",t(scale(t(mat[cor_atac_rna > .5,]))), cluster_columns = FALSE,  show_row_names = FALSE, column_split = c("ATAC","ATAC","ATAC","ATAC","ATAC","RNA","RNA","RNA","RNA","RNA"), column_title = "Positive correlation only")

draw(h1)
draw(h2)

