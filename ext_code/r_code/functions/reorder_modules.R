#wg_module <- smed_wg_module_
#order_criterion <- smed_ctypes$ctype
reorder_modules <- function(wg_module, order_criterion , min_kME = NULL, ordering_function = "upper quantile", thresh_sd = 1,...){
  module_col <- which(colnames(wg_module) == "module")
  
  if (!(is.null(min_kME))){
    print(paste0("using kME > ",min_kME))
    require(WGCNA)
    # Calculate eigengenes
    MEList <- moduleEigengenes(t(wg_module[,-module_col]), colors = wg_module[,module_col])
    MEs <- MEList$eigengenes
    datKME <- signedKME(t(wg_module[,-module_col]), MEs, outputColumnName = "MM_")
    
    filt_top <- 
      apply(
        datKME,
        1,
        function(x){
          if(any(x > min_kME)){
            res = TRUE
          } else {
            res = FALSE
          } 
          return(res)
        }
      )
    
    wg_module <- wg_module[filt_top,]
  }
  
  # pre-required function
  
  outliers <- function(x, thresh_sd){
    # get mean and Standard deviation
    m = mean(x)
    s = sd(x)
    
    # get threshold value for outliers
    Tmax = m+(thresh_sd*s)
    
    # find outlier
    ol <- which(x > Tmax)
    
    return(ol)
  }
  
  if(ordering_function == "upper quantile"){
    where_highest <- function(x){
      # get upper quantile for every cell type
      x_m <- sapply(x,quantile,0.75) 
      names(x_m) <- colnames(x)
      # which are the highest ones
      highest_filt <- outliers(x_m, thresh_sd = thresh_sd)
      # retrieve which ones and return
      y <- x_m[highest_filt]
      names(y) <- names(x_m)[highest_filt]
      y <- sort(y, decreasing = TRUE)
      return(y)
    }
  } else if(ordering_function == "median") {
    where_highest <- function(x){
      # get median for every cell type
      x_m <- sapply(x,median) 
      names(x_m) <- colnames(x)
      # which are the highest ones
      highest_filt <- outliers(x_m, thresh_sd = thresh_sd)
      # retrieve which ones and return
      y <- x_m[highest_filt]
      names(y) <- names(x_m)[highest_filt]
      y <- sort(y, decreasing = TRUE)
      return(y)
    }
  } else if(ordering_function == "mean") {
    where_highest <- function(x){
      # get mean for every cell type
      x_m <- sapply(x,mean) 
      names(x_m) <- colnames(x)
      # which are the highest ones
      highest_filt <- outliers(x_m, thresh_sd = thresh_sd)
      # retrieve which ones and return
      y <- x_m[highest_filt]
      names(y) <- names(x_m)[highest_filt]
      y <- sort(y, decreasing = TRUE)
      return(y)
    }
  } else {
    stop("Function not recognised.")
  }
  
  # separate table by modules and reorder
  l <- split(wg_module[,-module_col],wg_module$module)  
  l <- l[order(names(l))]
  
  nmodules <- length(l)
  
  l_hi <- lapply(
    l,
    where_highest
  )
  
  max_n <- max(sapply(l_hi,length))
  
  m_order <- matrix(0, nrow = length(l), ncol = max_n)
  rownames(m_order) <- names(l)
  
  for(i in 1:nmodules){
    r <- match(names(l_hi[[i]]),order_criterion)
    m_order[i,] <- 
      c(r,rep(0,max_n-length(r)))
  }
  
  m_order <-
    m_order[
      do.call(order,as.data.frame(m_order[,ncol(m_order):1])),
    ]
  
  res <- data.frame(
    neworder = 1:nmodules,
    module_wgcna = rownames(m_order),
    type = ifelse(apply(m_order, 1, function(x)length(which(x != 0))) == 1, "s", "m")
  )
  
  res$newname <- 
    unique(paste0(
      res$type,
      c(
        formatC(seq(1:length(which(res$type == "s"))), format = "d", width = 2, flag = "0"),
        formatC(seq(1:length(which(res$type == "m"))), format = "d", width = 2, flag = "0")
      )
    ))
  
  res$celltypes <-
    sapply(
      l_hi[match(res$module_wgcna,names(l_hi))],
      function(x){
        y = paste(names(x), collapse = ",")
        
        return(y)
      }
    )
  
  res$num_genes <-
    sapply(
      l[match(res$module_wgcna, names(l))],
      function(x)nrow(x)
    )
  
  
  return(res)
  }

