deseq_sc <- 
  function(
    m,d,
    filter_by = "padj", p_threshold = 0.05, lfc_threshold = 0, contrast_info,
    plot_results = TRUE, cell, min_passing_samples = 3, min_counts_per_sample = 3,
    keep_dubious = FALSE, ...
  ) {
    
    # Setup
    if( !(filter_by %in% c("padj","pvalue"))) 
      stop("E: filtering column specified is not 'pvalue' or 'padj'.")
    
    require(DESeq2)
    require(dplyr)
    require(ggplot2)
    require(ComplexHeatmap)
    
    #' Params:
    #' @param x: a gene expression matrix with samples in columns
    quantnorm <- function(x){
      require(preprocessCore)
      colnames_x <- colnames(x)
      rownames_x <- rownames(x)
      x <- as.matrix(x)
      x_q <- normalize.quantiles(x) # figure out why we get zeroes and neg numbers here
      colnames(x_q) <- colnames_x
      rownames(x_q) <- rownames_x
      return(x_q)
    }
    
    #' Params:
    #' @param x: a gene expression matrix with samples in columns.
    #' @param tbl: a sample table with colnames(x) in "sample" column, 
    #' condition info (which will be the col names of output exp matrix) in "condition",
    #' and replicate info in "rep" conditions.
    sum_by_repl <- function(x,tbl){
      # create matrix where to put counts
      m <- matrix(0, nrow = nrow(x),ncol = length(unique(tbl$condition)))
      rownames(m) <- rownames(x)
      colnames(m) <- unique(tbl$condition)
      
      for (i in colnames(m)){
        filt <- which(colnames(x) %in% tbl$sample[tbl$condition == i])
        x_i <- x[,filt]
        
        if(length(filt) > 1){
          m[,i] <- rowSums(x_i)
        } else {
          m[,i] <- x_i
        }
      }
      
      return(m)
    }
    
    # Setup Variables
    conds <- contrast_info[c(3,2)]
    
    cond_col <- which( colnames(d) == contrast_info[1] )
    
    d <- d[d[[cond_col]] %in% conds,]
    d[[cond_col]] <- droplevels(d[[cond_col]])
    d <- d[,sapply(d,function(x)length(unique(x))) != 1]
    
    m <- m[,colnames(m) %in% d[,1]]
    
    # DESeq object
    dds <- DESeqDataSetFromMatrix(
      countData = m,
      colData = d,
      design = formula(paste("~","replicate","+",contrast_info[1] ) ), # design = ~replicate + condition, #formula(paste("~",x)))
      ...
    )
    
    print("filtering")
    smallestGroupSize <- min_passing_samples
    keep <- rowSums(counts(dds) >= min_counts_per_sample) >= smallestGroupSize
    # (we could expand the list of genes we are not recapturing using the following)
    if(keep_dubious == TRUE) {
      d2 <- d
      d2$condition == d[[cond_col]]
      keep2 <- apply(sum_by_repl(counts(dds),d2) >= min_counts_per_sample,1,function(x) any(x == TRUE))
      keep[which(!(names(keep[keep == TRUE]) %in% names(keep2[keep2==TRUE])))] <- TRUE
    }
    
    dds <- dds[keep,]
    
    print("DGE")
    dds <- DESeq(dds, ...)
    
    # Results
    print("results")
    res <- results(object = dds, contrast = contrast_info)
    r <- as.data.frame(res)
    
    # Retrieve diff genes
    filter_by_column <- which(colnames(r) == filter_by)
    signif <- rownames(r)[r[,filter_by_column] < p_threshold]
    signif <- signif[complete.cases(signif)]
    
    # Special cases if 1 or 0 diffgenes found
    if(length(signif) < 1 ){
      if (length(signif) == 1) {
        res <- list(
          dds <- dds,
          res = r,
          diffreg_genes = signif
        )
        print("only 1 diff gene found")
        return(res)
      } else {
        res <- list(
          dds = dds,
          res = r,
          diffgenes = NA
        )
        print("No diff genes found")
        return(res)
      }
    }
    
    # Build results list
    res <- list(
      dds = dds,
      res = r,
      diffgenes = signif
    )
    
    # If plotting has been requested
    if(plot_results == TRUE){
      # Create volcano plot
      print("volcano")
      v <-
        r %>% 
        ggplot(
          aes(x = log2FoldChange, y = -log(r[[filter_by]]))
        ) + 
        geom_point() + theme_minimal()
      
      # Create transformed data matrix for heatmap
      print("heatmap qnorm")
      m_hm <-  as.matrix(t(scale(t(quantnorm(m)[rownames(m) %in% signif,]))))
      
      # Create heatmap annotations
      print("hm annot")
      ha_sampletable <- HeatmapAnnotation(
        df = d[,-c(1)], 
        col = dynamic_colors_annotation(df = d[,-c(1)], rand.seed = 4)
      )
      
      # Create heatmap
      print("hm")
      hm <- Heatmap(
        name = "expression",
        m_hm,
        clustering_method_rows = "ward.D2",
        bottom_annotation = ha_sampletable,
        column_title = cell
      )
      
      # Create pca
      print("vsd")
      vsd <-  varianceStabilizingTransformation(dds, blind = FALSE) # already filtered
      print("pca")
      pcaData <- plotPCA(vsd, intgroup=c(contrast_info[1], "replicate"), returnData=TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      pc <- 
        ggplot(pcaData, aes(PC1, PC2, color=.data[[contrast_info[1]]], shape=replicate)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed()
      plots <- list(
        volcano = v,
        heatmap = hm,
        pca = pc
      ) 
      res <- append(res,plots)
    }
    
    # Finish and return results
    print(paste0("Found ",length(signif)," diff genes."))
    return(res)
  }