###
#Functions are adopted and modified from @alberto



#WGCNA_prep_pseudobulk
wgcna_prep=function(obj,group,ctype_info_path,mxtype,hclust.md){
  
  obj= subset(obj, orig.ident == group)
  dir.create(paste0("./figures/",group,"/"))
  dir.create(paste0("./figures/",group,"/WGCNA"))
  dir.create(paste0("./figures/",group,"/WGCNA/wgcna_exploration_modules"))
  dir.create(paste0("./figures/",group,"/WGCNA/wgcna_exploration_modules/featplots_modules_deseq_CW"))
  dir.create(paste0("./outputs/",group,"/"))
  dir.create(paste0("./data/",group,"/"))
  dir.create(paste0("./data/",group,"/rda/"))
  dir.create(paste0("./data/",group,"/rds/"))
  
  ## Getting cluster information for pseudobulk
  
  #scRNA-seq:
  DefaultAssay(obj)<-"RNA"
  ss_clusters <- 
    setNames(
      obj$merged_sub.anno_type,
      colnames(obj)
    )
  
  ss_pseudobulk <- 
    pseudobulk(
      x = obj@assays$RNA$counts,
      ident = ss_clusters
    )
  
  
  #{r eval = FALSE}
  ss_psbulk_ncells <- 
    pseudobulk_ncells(
      x = obj@assays$RNA@counts,
      identities = Idents(obj),
      min_counts = 1
    )
  
  cluster_size = 
    c(table(obj$merged_sub.anno_type))
  
  
 
  #reorder column
  coltb=data.frame("id"=obj$merged_sub.anno_type,
                   "col"=obj$ctype_col)
  coltb=arrange(coltb,id)
  coltb=unique(coltb)
  #barclaobj=sapply(str_split(colnames(ss_pseudobulk),"_"), `[`, 1) 
  
  cols=translate_ids(coltb$id,dict = coltb[,c(1,2)])
  
  a=data.frame("n.gene"=rev(sapply(as.data.frame(ss_pseudobulk), function(x){length(which(x > 5))})))
  a$col=rev(cols)
  a=a[c(str_sort(colnames(ss_pseudobulk),numeric = T)),]
  ag=a$n.gene
  names(ag)= rownames(a)
  
  pdf(paste0("./figures/",group,"/supp_1_1_G.pdf"), wi = 6, he = 6)
  par(mar = c(5,10,2,6)+.1)
  barplot(
    rev(ag),  col = alpha(rev(a$col),0.8),
    border = darken(rev(a$col), 0.6),
    las = 1, horiz = TRUE, 
    cex.names = .3,
    xlim = c(0,15000),
    xlab = "No. genes quantified on each cluster (>5 counts)"
  )
  dev.off()
  
  ## Save Data
  
  
  write.table(
    ss_pseudobulk,
    file = paste0("./data/",group,"/ss_pseudobulk.tsv"),
    sep = "\t", dec = ".", 
    row.names = TRUE, quote = FALSE
  )
  
  write.table(
    ss_psbulk_ncells,
    file = paste0("./data/",group,"/ss_pseudobulk_ncells.tsv"),
    sep = "\t", dec = ".", 
    row.names = TRUE, quote = FALSE
  )
  
  write.table(
    cluster_size,
    file = paste0("./data/",group,"/ss_clustersize.tsv"),
    sep = "\t", dec = ".", 
    row.names = TRUE, quote = FALSE
  )
  

  # 09 ----------------------------------------------------------------------
  library(vroom)
  library(reshape2)
  library(ComplexHeatmap)
  library(circlize)
  library(viridis)
  library(colorspace)
  library(RColorBrewer)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(gplots)
  library(data.table)
  library(DESeq2)
  library(Matrix)
  
  
  ## Loading cell type information
  
  
  danio_ctypes <- read.csv(ctype_info_path,row.names = 1)
  danio_ctypes$ctype <- factor(danio_ctypes$ctype,levels = unique(danio_ctypes$ctype))
  danio_ctypes=danio_ctypes[c(rownames(obj[["merged_sub.anno_type"]])),]
  danio_ctypes$merged_sub.anno_type=obj$merged_sub.anno_type
  danio_ctypes$ctype <- factor(danio_ctypes$merged_sub.anno_type,levels = unique(danio_ctypes$merged_sub.anno_type))
  
  ## Loading counts
  
  
  danio_counts <- ss_pseudobulk
  danio_counts <- danio_counts[rowSums(danio_counts) > 10,]
  
  
  #And if we check the distribution of counts per cluster:
  
  
  psbulk_ncells <- ss_psbulk_ncells
  
  
  #We will calculate a cell weight matrix to weigh the expression values from the pseudobulk count matrix. For this it will calculate how many cells (in percentage) are expressing a given gene in a given cluster, in relation to how many cells (in percentage) are expressing that gene in the rest of clusters. We set the minimum number of total counts as 30 and the minimum cells to take into account as 1.
  
  
  min_counts <- 30
  min_cells <- 3
  cluster_size = cluster_size
  
  
  #Now we create the weight values matrix:
  
  
  danio_psbulk_cellweights <- 
    get_cellweight_matrix(
      x = as.matrix(danio_counts),
      y = psbulk_ncells,
      C = cluster_size,
      min_counts = min_counts,
      min_cells = min_cells
    )
  
  
  #A quick look at this matrix:
  
  
  danio_psbulk_cellweights[1:5,1:5]
  
  dim(danio_psbulk_cellweights)
  
  
  #We subset the counts matrix to keep the same genes that were retrieved in the cell weight values matrix
  colnames(danio_counts)=gsub("\\/",".",colnames(danio_counts))
  danio_ctypes$ctype=gsub("\\/",".",danio_ctypes$ctype)
  
  
  m_ <- as.matrix(danio_counts)
  m_ <- m_[rowSums(m_) >= min_counts,]
  m_ <- m_[rownames(m_) %in% rownames(danio_psbulk_cellweights),]
  
  
  #And we normalise using DESeq2
  
  
  m_dds <- DESeqDataSetFromMatrix(
    countData = m_,
    colData = data.frame(condition = colnames(m_)),
    design = ~ condition)
  m_dds <- estimateSizeFactors(m_dds)
  danio_counts_norm <- counts(m_dds, normalized=TRUE)
  
  
  #The final matrix is the log-transformed of these DESeq2-normalised values, multiplied by the cell weights.
  danio_counts_norm_cw <- (log1p(danio_counts_norm) * danio_psbulk_cellweights)
  
  #And here the boxplots:
  
  #{r, fig.width = 8, fig.height = 6}
  # Boxplots of count normalisation
  
  tiff(paste0("./figures/",group,"/psudo_bulk.tiff"),
       width = 40,height = 30,units = "cm", res = 300,compression = "lzw")
  par(mfrow = c(2,2))
  
  boxplot(
    log1p(m_),
    las = 2,
    cex.axis=0.5,
    cex = 0.5,
    col = danio_ctypes$ctype_col,
    border = darken(danio_ctypes$ctype_col),
    pch = 16,
    outcol=rgb(0.1,0.1,0.1,0.1),
    main = "raw"
  )
  
  boxplot(
    log1p(danio_counts_norm),
    las = 2,
    cex.axis=0.5,
    cex = 0.5,
    col = danio_ctypes$ctype_col,
    border = darken(danio_ctypes$ctype_col),
    pch = 16,
    outcol=rgb(0.1,0.1,0.1,0.1),
    main = "raw, post-deseq2"
  )
  
  boxplot(
    log1p(m_)*danio_psbulk_cellweights,
    las = 2,
    cex.axis=0.5,
    cex = 0.5,
    col = danio_ctypes$ctype_col,
    border = darken(danio_ctypes$ctype_col),
    pch = 16,
    outcol=rgb(0.1,0.1,0.1,0.1),
    main = "ncell-informed"
  )
  
  boxplot(
    danio_counts_norm_cw,
    las = 2,
    cex.axis=0.5,
    cex = 0.5,
    col = danio_ctypes$ctype_col,
    border = darken(danio_ctypes$ctype_col),
    pch = 16,
    outcol=rgb(0.1,0.1,0.1,0.1),
    main = "ncell-informed, post-deseq2"
  )
  
  par(mfrow = c(1,1))
  dev.off()
  
  ## Assigning every gene to a cell type
  
  #This is a quick way to assign a gene to every cell type by simply pinning what is the cell type with the highest expression of a given gene.
  danio_counts_norm_cw=get(mxtype)
  
  danio_genecolor <- data.frame(
    id = rownames(danio_counts_norm_cw),
    genecolor = apply(
      danio_counts_norm_cw,
      1,
      function(x) {
        a <- which( x == max(x) )  
        b <- names(x[a])
        return(b) # assign in the table
      }
    )
  )
  
  danio_genecolor$genecolor <- translate_ids(danio_genecolor$genecolor, dict = danio_ctypes[,c(1,3)])
  
  head(danio_genecolor)
  
  
  ## Defining broad cell types: co-occurrence of cell cluster similarity
  
  #We can use a survival clustering approach (Levy et al., 2021) to infer what are the most similar cell type clusters.
  #We will apply a soft threshold for genes with CV higher than 0.25.
  
  
  # comparABle function tidyup from source
  source("./ext_code/comparABle-main/comparABle-main/code/functions/tidyup_functions.R")
  
  danio_cpm_cooc <-
    tidyup(
      danio_counts_norm_cw,#[rownames(danio_cpm) %in% danio_hvgs,],
      highlyvariable = TRUE #FALSE # TRUE if not using the subset of danio_hvgs
    )
  
  # set fixed seed
  set.seed(4343)
  h <- c(0.75,0.9)
  clustering_algorithm <- "hclust"
  clustering_method <- hclust.md 
  cor_method <- "pearson"
  p <- 0.01
  danio_cpm_vargenes = rownames(danio_cpm_cooc)
  
  # Levy et al 2021 'treeFromEnsembleClustering' from source
  source("./ext_code/r_code/functions/treeFromEnsembleClustering.R")
  cooc <- treeFromEnsembleClustering(
    x=danio_cpm_cooc, p=p, h=h,  n = 10000, vargenes = danio_cpm_vargenes, bootstrap=FALSE,
    clustering_algorithm=clustering_algorithm, clustering_method=clustering_method, 
    cor_method=cor_method
  )
  
  
  #The resulting heatmap of similarity:
  
  #{r fig.width=7.5, fig.height=6.5, echo = FALSE}
  ctypes_rowAnno <-
    rowAnnotation(
      cluster = rownames(cooc$cooccurrence),
      col = list( cluster = setNames(danio_ctypes$ctype_col ,danio_ctypes$ctype) ),
      show_legend = F, show_annotation_name = F
    )
  
  
  clu_ha = HeatmapAnnotation(
    name = "cell types",
    cluster = factor(colnames(cooc$cooccurrence), levels = unique(danio_ctypes$ctype)),
    col = list(cluster = setNames(danio_ctypes$ctype_col ,danio_ctypes$ctype)),
    show_legend = F, show_annotation_name = F
  )
  
  
  #mts=cooc$cooccurrence
  #mts=mts[str_sort(rownames(mts),numeric = T),str_sort(colnames(mts),numeric = T)]
  danio_cooc_hm <- Heatmap(
    name="co-occurence",
    cooc$cooccurrence,
    col = colorRamp2(
      c(seq(min(cooc$cooccurrence),
            max(cooc$cooccurrence),
            length=9
      )
      ),
      colors=c( 
        c("#FFFFEA","#ffffe5","#fff7bc","#fee391","#fec44f","#fe9929","#ec7014","#cc4c02","#990000")
      )
    ),
    cluster_rows = as.hclust(cooc$tree),
    cluster_columns = as.hclust(cooc$tree),
    left_annotation = ctypes_rowAnno,
    top_annotation = clu_ha,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8)
  )
  
  draw(danio_cooc_hm)
  
  tiff(paste0("./figures/",group,"/WGCNA/cooccurrence.tiff"),
       width = 45,height = 40,units = "cm", res = 300,compression = "lzw")
  draw(danio_cooc_hm)
  
  dev.off()
  
  ## Save the data
  
  
  save(
    danio_ctypes,
    danio_counts,
    danio_psbulk_cellweights,
    danio_counts_norm,
    danio_counts_norm_cw,
    danio_genecolor,
    file = paste0("./data/",group,"/rda/danio_counts.rda")
  )
  
  
}

#WGCNA running
#sp=sps[i],hclust.md=hclust.method,min.md.size=30,group=libr[i]
run_WGCNA = function(obj,sp,hclust.md,min.md.size,group){
  load(paste0("./data/",group,"/rda/danio_counts.rda"))
  
  softPower  <- sp
  
  
  ## The Adjacency Matrix
  
  "WGCNA's criterion, for two genes to be adjacent, is that they show similar levels of co-regulation. Co-regulation is defined high values of signed Pearson Correlation (that is, they show high levels of correlation or anti-correlation).

The output is a matrix of n genes x n genes showcasing the level of corregulation between genes.
"
  #{r adjacency, eval = TRUE, echo = TRUE}
  adjacency  <-  adjacency(datExpr, power = softPower)
  dim(adjacency)
  
  
  ## The Topological Overlapped Matrix (TOM)
  '
But WGCNA does not leave it there. This method takes into account the amount of shared neighbourhoods between pairs of genes to strengthen or weaken the level of association between the two.

This is ran on top of the adjacency matrix and can take a LOT of time. For practicality, we have included the TOM we generated as an .rda object that we load to keep up with the analysis.

As before, the TOM matrix is a n genes x n genes matrix.
'
  #{r TOM}
  TOM <- TOMsimilarity(adjacency)
  dimnames(TOM) <- dimnames(adjacency)
  dim(TOM)
  
  '
The TOM can be used to generate graph objects that can be analysed using igraph. We will do that later.

Higher values of TOM indicate higher level of association. This metric can be transformed to depict (dis)similarity between the genes, which gets us closer to the clustering steps.
'
  #{r dissTOM}
  dissTOM <- 1 - TOM
  
  
  ## Clustering and module detection
  '
WGCNA uses hierarchical clustering of the TOM-derived similarity to determine module membership.
'
  #{r cluster_genes_by_dissTOM}
  geneTree = hclust(as.dist(dissTOM), method = hclust.md)
  
  '
We set a minimum module size of thirty genes, and ask to cut the tree. See `?cutreeDynamic` for a larger explanation.

Importantly: The output value is a vector of numerical labels giving assignment of objects to modules. Unassigned objects are labeled 0, the largest module has label 1, next largest 2 etc.
'
  #{r assign_genes}
  minModuleSize <- min.md.size
  
  dynamicMods <- cutreeDynamic(
    dendro = geneTree,
    distM = dissTOM,
    deepSplit = 3,
    pamRespectsDendro = FALSE,
    minClusterSize = minModuleSize
  )
  
  
  #The output vector, can be transformed from discrete numeric values to colors for a more qualitative, but equally blindfolded naming system. Because the order of the contents of the vector matches the order of the genes in the datExpr matrix and the TOM matrix, we can transfer the gene names to this vector to retrieve an association gene-module.
  
  #{r rename_modules_as_colors}
  moduleColors <- labels2colors(dynamicMods)
  names(moduleColors) <- colnames(datExpr)
  
  
  #And even more, we can transform this information in a more human-friendly (and dplyr-friendly) format that we will use extensively in our downstream analysis:
  
  #{r data_frame_danio_modules}
  danio_id_module_wgcna <- data.frame(
    id = colnames(datExpr),
    module = moduleColors
  )
  
  
  #We can visualise the looks of our modules with `plotDendroAndColors`.
  
  #{r plot_dendro_and_colors}
  plotDendroAndColors(
    geneTree, moduleColors,
    c("Dynamic Tree Cut"), dendroLabels = FALSE,
    hang = 0.03, addGuide = TRUE, guideHang = 0.05
  )
  
  #{r}
  length(unique(moduleColors))
  
  tiff(paste0("./figures/",group,"/WGCNA/dendrogram.tiff"),
       width = 20,height = 6,units = "cm", res = 300,compression = "lzw")
  plotDendroAndColors(
    geneTree, moduleColors,
    c("Dynamic Tree Cut"), dendroLabels = FALSE,
    hang = 0.03, addGuide = TRUE, guideHang = 0.05
  )
  dev.off()
  
  
  '
Based on the early splitting and length of branches, from this plot we can observe that gene modules are defined very discretely. We were able to observe a similar trend when browsing the heatmaps of the transcription factors in markdown #02.

This can change for every dataset based on the species, quality of the data, and overall conditions of the experiment.

WGCNA allows for extra steps to merge together modules that are too similar that were perhaps accidentally split when cutting the tree. This was not done for Schmidtea as gene modules proved to be very modular.

As said before, refer to the official documentation for a more detailed depiction of the step-by-step-analysis.
'
  ## Renaming and reordering the modules
  
  #{r}
  # Prep
  danio_wg_module <- 
    merge(
      t(datExpr),
      danio_id_module_wgcna,
      by.x = 0, by.y = 1,
      all.X = TRUE
    )
  rownames(danio_wg_module) <- danio_wg_module[,1]
  danio_wg_module[,1] <- NULL
  
  
  #{r}
  # reordered modules
  danio_modules_table <- reorder_modules(danio_wg_module, 
                                        order_criterion =danio_ctypes$ctype,ordering_function = "median", 
                                        thresh_sd = 1.5)
  
  danio_modules_table$newname <- factor(danio_modules_table$newname,
                                       levels=unique(danio_modules_table$newname))
  
  
  danio_modules_table$newcolor <- 
    translate_ids(
      x = danio_modules_table$module_wgcna,
      dict = 
        unique(data.frame(
          module_wgcna = danio_modules_table$module_wgcna,
          colour = sapply(as.data.frame(rgb2hsv(col2rgb(danio_modules_table$module_wgcna))),pastelise_hsv,n=0.6)
        ))
    )
  
  
  #{r}
  danio_id_module <-
    data.frame(
      id = danio_id_module_wgcna$id,
      module = translate_ids(x=danio_id_module_wgcna$module,dict = danio_modules_table[,c(2,4)])
    )
  head(danio_id_module)
  
  
  hm_bp <-
    t(apply(
      aggregate(danio_counts_norm_cw[danio_id_module$id,],by = list(danio_id_module$module),FUN = mean)[,-1],
      1,
      function(x)x/sum(x)
    ))
  rownames(hm_bp) <- levels(danio_modules_table$newname)
  
  danio_modules_table$cell_color <-
    translate_ids(
      x = as.character(danio_modules_table$newname),
      dict = data.frame(
        module = rownames(hm_bp),
        color = translate_ids(apply(hm_bp,1,function(x){y=names(x[x==max(x)]); return(y)}), danio_ctypes[,c(1,3)])
      )
    )
  
  head(danio_modules_table)
  
  
  
  #{r}
  danio_wg_module$module <- 
    factor(
      translate_ids(danio_wg_module$module, dict = danio_modules_table[,c(2,4)]),
      levels = unique(danio_modules_table$newname)
    )
  
  danio_wg_module  <- danio_wg_module[order(danio_wg_module$module),]
  
  
  
  
  #{r}
  s_m_distrs = 
    reshape2::melt(danio_wg_module) %>%
    mutate(celltype = factor(variable,levels = unique(danio_ctypes$ctype))) %>% 
    group_by(module,celltype) %>% dplyr::summarise(up_q = quantile(value,.75)) %>% 
    mutate(norm_up_q = relativise(up_q))
  
  p_s_m_distrs = 
    s_m_distrs %>%
    ggplot(aes(x=norm_up_q,y=module,col=celltype))+
    geom_boxplot(outliers = FALSE)+
    geom_jitter(position = "dodge")+
    scale_color_manual(values = alpha(danio_ctypes$ctype_col,.75))+
    scale_y_discrete(limits = rev(levels(danio_wg_module$module)))+
    ylab("Module")+
    xlab("Normalised upper quantile of expression, on each cell type")+
    theme_classic()+
    theme(legend.position="none")+
    ggtitle("Gene expression dynamics of each module")
  
  print(p_s_m_distrs)
  
  
  pdf(paste0("./figures/",group,"/WGCNA/s_m_distributions_plot.pdf"), he = 8, wi = 3)
  print(p_s_m_distrs)
  dev.off()
  
  ## Barplots of all the modules
  
  #{r}
  # Barplots of all the gene modules
  modulecolumn <- which(!sapply(danio_wg_module,is.numeric))
  pdf(
    file = paste0("./figures/",group,"/WGCNA/danio_wgcna_module_boxplots.pdf"),
    width = 10, height = 5
  )
  par(
    mar=c(12,4,4,2)+0.1,
    xpd = TRUE
  )
  
  for (i in sort(unique(danio_wg_module$module))) {
    module_i <- rownames(danio_wg_module)[danio_wg_module$module == i]
    module_i_counts <- danio_wg_module[rownames(danio_wg_module) %in% module_i,-modulecolumn]
    
    boxplot(
      log1p(module_i_counts),
      col = danio_ctypes$ctype_col,
      border = darken(danio_ctypes$ctype_col,0.6),
      las=2,
      cex = 0.5,
      cex.axis=0.5,
      pch = 16,
      ylab = "log1p(counts)",
      outcol=rgb(0.1,0.1,0.1,0.1),
      main = paste0("Module ",i,"; ngenes: ", length(module_i))
    )
  }
  dev.off()
  ## Feature plots of all the modules
  
  #{r}

  
  #{r, message=FALSE}
  # Featureplots of all gene modules
  featplot_genes_from_modules(
    module_list = danio_modules_table$newname,
    wg_module = danio_wg_module,
    scdata = obj,
    output_root_path = paste0("./figures/",group,"/WGCNA/wgcna_exploration_modules/featplots_modules_deseq_CW"),
    num_genes=50
  )
  
  
  ## Plotting the network
  
  #{r tom_network_plot, echo = TRUE}
  restGenes <- (moduleColors != "grey")
  diss <- 1-TOMsimilarityFromExpr( datExpr[, restGenes], power = softPower )
  hier1 <- hclust(as.dist(diss), method=hclust.md )
  diag(diss) = NA
  
  pdf(
    file = paste0("./figures/",group,"/WGCNA/danio_wgcna_NetworkTOMplot.pdf"),
    width = 10,
    height = 10
  )
  TOMplot(1-diss^4, hier1, as.character(moduleColors[restGenes]), # 1-X to change color
          main = "danio TOM heatmap plot, module genes" )
  dev.off()
  
  
  #And the ouptut:
  
  #(insert output here)
  
  ## Module eigengenes and connectivity
  
  #{r}
  # Calculate eigengenes
  MEList <- moduleEigengenes(datExpr, colors = moduleColors,excludeGrey = F)
  MEs <- MEList$eigengenes
  MEs <- MEs[
    ,
    match(
      danio_modules_table$newname,
      translate_ids(gsub("ME","",colnames(MEs)),dict = danio_modules_table[,c(2,4)])
    )]
  
  colnames(MEs) <- 
    paste0(
      "ME",
      as.character(
        translate_ids(
          gsub(".1","",gsub("ME","",colnames(MEs))),
          dict = danio_modules_table[,c(2,4)]
        )
      )
    )
  
  datKME <- signedKME(datExpr, MEs, outputColumnName = "")
  
  min_kme <- 0.9
  
  filt_top <- 
    apply(
      datKME, 1,
      function(x){
        if(any(x > min_kme)){
          res = TRUE
        } else {
          res = FALSE
        } 
        return(res)
      }
    )
  
  danio_wg_module_top <- danio_wg_module[filt_top,]
  
  
  ## Visualisation
  
  #A tidier version of our expression data, organised by module membership:
  
  #{r}
  set.seed(4343)
  danio_wg_module_viz <- danio_wg_module_top %>% group_by(module) %>% slice_sample(n=30)
  danio_wg_module_viz <- danio_wg_module_viz[complete.cases(danio_wg_module_viz),]
  
  #{r, fig.height=12, fig.width=4, message=FALSE, warning=FALSE}
  modulecolumn <- which(!sapply(danio_wg_module,is.numeric))
  danio_wg_module_viz=danio_wg_module_viz[,c(str_sort(colnames(danio_wg_module_viz),numeric = TRUE))]
  
  str_sort(unique(danio_ctypes$ctype),numeric = TRUE)
  
  setcol=danio_ctypes[,c("ctype_col","ctype")]
  setcol=unique(setcol)
  setcol=dplyr::arrange(setcol,ctype)
  clu_ha = HeatmapAnnotation(
    name = "cell types",annotation_height = 1,
    cluster = factor(setcol$ctype, levels = setcol$ctype),
    col = list(cluster = setNames(setcol$ctype_col,setcol$ctype)),
    show_legend = FALSE
  )
  
  mt=as.matrix(danio_wg_module_viz[,-modulecolumn])
  mt=mt[,c(setcol$ctype)]
  danio_wg_hm <- ComplexHeatmap::Heatmap(
    name = "z-score",
    mt, #danio_wg_module_viz[,-modulecolumn], 
    cluster_rows= F,
    show_row_names = F,
    show_row_dend = F,
    cluster_columns = F,
    show_column_names = TRUE,
    column_names_side = "bottom",
    row_split = danio_wg_module_viz$module,
    row_title_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 6),
    row_title_side = "left",
    row_title_rot = 0,
    top_annotation = clu_ha,
    bottom_annotation = clu_ha,
    use_raster = FALSE,
    col = colorRampPalette(c("#f1f5ff","#b1b8c4","#2c2b46"))(10)
  )
  draw(danio_wg_hm)
  
  
  #{r}
  pdf(paste0("./figures/",group,"/WGCNA/wgcna_exploration_modules/danio_wgcna_module_exploration_heatmap.pdf"), height = 15, width = 5)
  draw(danio_wg_hm)
  dev.off()
  
  mt=as.matrix(danio_wg_module_viz[,-modulecolumn])
  
  danio_wg_hm <- ComplexHeatmap::Heatmap(
    name = "z-score",
    mt, #danio_wg_module_viz[,-modulecolumn], 
    cluster_rows= F,
    show_row_names = F,
    show_row_dend = F,
    cluster_columns = F,
    show_column_names = TRUE,
    column_names_side = "bottom",
    row_split = danio_wg_module_viz$module,
    row_title_gp = gpar(fontsize = 10),
    #column_names_gp = gpar(fontsize = 6),
    row_title_side = "left",
    row_title_rot = 0,
    top_annotation = clu_ha,
    bottom_annotation = clu_ha,
    use_raster = FALSE,
    col = colorRampPalette(c("#f1f5ff","#b1b8c4","#2c2b46"))(10)
  )
  draw(danio_wg_hm)
  
  
  #{r}
  pdf(paste0("./figures/",group,"/WGCNA/wgcna_exploration_modules/danio_wgcna_module_exploration_heatmap2.pdf"), height = 15, width = 5)
  draw(danio_wg_hm)
  dev.off()
  
  #And the result of plotting this as a heatmap, where rows correspond to genes and columns correspond to cell clusters. Color intensity indicates z-scored expression at a given cell cluster. Several transcription factors of interest have been highlighted.
  
  #{r echo = FALSE, warning = FALSE, message = FALSE, fig.width = 6, fig.height = 15}
  ngenes_per_module <- 30
  
  set.seed(4343)
  danio_wg_module_viz <- danio_wg_module %>% group_by(module) %>% slice_sample(n=ngenes_per_module)
  modulecolumn2 <- which(!sapply(danio_wg_module_viz,is.numeric))
  
  setcol=danio_ctypes[,c("ctype_col","ctype")]
  setcol=unique(setcol)
  setcol=dplyr::arrange(setcol,ctype)
  
  wg_ha = HeatmapAnnotation(
    name = "cell types",
    cluster = factor(setcol$ctype, levels=setcol$ctype),
    col = list(
      cluster = setNames(
        setcol$ctype_col, # this was here before just in case needed again [match(colnames(danio_wg_module),danio_ctypes$ctype)]
        setcol$ctype) # same
    ),
    show_legend = FALSE
  )
  
  txt = lapply(danio_modules_table$newname,FUN=function(x){x})
  names(txt) <- unique(levels(danio_id_module$module))
  
  moduleSizes <- setNames(danio_modules_table$num_genes,levels(danio_modules_table$newname))
  
  wg_mod <-
    HeatmapAnnotation(
      stack = 
        anno_barplot(
          hm_bp[rep(1:length(levels(danio_modules_table$newname)),each=ngenes_per_module),setcol$ctype],
          gp = gpar(col = setcol$ctype_col), # this was here [match(colnames(danio_wg_module)[1:36],danio_ctypes$ctype)][-37]
          border = FALSE
        ),
      log10_ngen = 
        anno_barplot(
          log10(rep(moduleSizes,each=ngenes_per_module)),
          gp = gpar(col = "#444444"),#rep(danio_modules_table$general_color,each=30),xlim = c(50,1000)),
          border = FALSE
        ),
      gap = unit(10, "points"),
      which = "row"
    )
  wg_ha@anno_list$cluster@label <- NULL
  
  mt= danio_wg_module_viz[,-modulecolumn2]+2
  mt=mt[,setcol$ctype]
  danio_wgcna_hm <-
    Heatmap(
      mt,
      name = "expression",
      cluster_rows= F,
      show_row_names = F,
      show_row_dend = F,
      cluster_columns = F,
      show_column_names = TRUE,
      column_names_side = "bottom",
      row_split = danio_wg_module_viz$module,
      row_title_gp = gpar(fontsize = 10),
      column_names_gp = gpar(fontsize = 6),
      row_title_side = "left",
      row_title_rot = 0,border_gp = gpar(color = "grey"),#column_split = 
      col = colorRampPalette(c("#f1f5ff","#b1b8c4","#2c2b46"))(10),
      top_annotation=wg_ha,
      bottom_annotation = wg_ha,
      right_annotation = wg_mod,
      heatmap_legend_param = gpar(nrow = 2)
    )
  
  draw(danio_wgcna_hm)
  
  
  #{r}
  pdf(paste0("./figures/",group,"/WGCNA/danio_wgcna_heatmap.pdf"), width = 8, heigh = 12)
  draw(danio_wgcna_hm)
  dev.off()
  
  
  #{r save_id_module}
  write.xlsx(
    danio_id_module,
    file = paste0(
      "./outputs/",group,"/",
      "danio_wgcna_id_module.xlsx"
    ),
    sheetName = "schmidtea_wgcna_id_module",
    col.names = TRUE,
    row.names = FALSE,
    showNA = TRUE
  )
  
  
  #We will store a second TOM, pruned in such a way, to generate graphs in the upcoming analyses.
  
  #{r save_TOM_graph_analysis}
  TOM_2 <- TOM[
    rownames(TOM) %in% danio_id_module$id,
    colnames(TOM) %in% danio_id_module$id
  ]
  saveRDS(TOM_2,paste0("data/",group,"/rda/wgcna_TOM_matrix.rds"))
  
  
  #{r save_all_wgcna}
  save(
    danio_modules_table,
    danio_id_module,
    datExpr,
    danio_wg_module,
    datKME,
    #danio_wg_GO_all,
    geneTree,
    MEList,
    MEs,
    #sft,
    danio_wg_module_viz,
    #danio_wg_list,
    adjacency,
    TOM,
    dissTOM,
    TOM_2,
    hm_bp,
    file = paste0("data/",group,"/rda/danio_wgcna_all.rda")
  )
  
  
  #{r save_subset_of_wgcna_for_graph_analysis}
  save(
    danio_modules_table,
    danio_id_module,
    datExpr,
    danio_wg_module,
    datKME,
    MEs,
    #danio_wg_GO_all,
    hm_bp,
    file = paste0("data/",group,"/rda/danio_wgcna.rda")
  )
  
  xlsx::write.xlsx(
    danio_modules_table,
    file = paste0("outputs/",group,"/ss_modules_table.xlsx"),
    sheetName = "module information",
    col.names = TRUE,
    row.names = FALSE
  )
  
  
  
  #{r}
  danio_wg_hm_complete <- Heatmap(
    name = "z-score",
    danio_wg_module[,-modulecolumn],
    cluster_rows= F,
    show_row_names = F,
    show_row_dend = F,
    cluster_columns = F,
    show_column_names = TRUE,
    column_names_side = "bottom",
    row_split = danio_wg_module$module,
    row_title_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 6),
    row_title_side = "left",
    row_title_rot = 0,
    top_annotation = clu_ha,
    bottom_annotation = clu_ha,
    use_raster = FALSE,
    col = colorRampPalette(c("#f1f5ff","#b1b8c4","#2c2b46"))(10)
  )
  
  pdf(paste0("./figures/",group,"/WGCNA/danio_wgcna_heatmap_complete.pdf"), wi = 3, he = 15)
  draw(danio_wg_hm_complete)
  dev.off()
  
  pdf(paste0("./figures/",group,"/WGCNA/danio_wgcna_dendro_and_colors_and_legend.pdf"), width = 10, height = 6)
  plotDendroAndColors(
    geneTree, translate_ids(danio_id_module$module,danio_modules_table[,c(4,7)]),
    c("Dynamic Tree Cut"), dendroLabels = FALSE,
    hang = 0.03, addGuide = TRUE, guideHang = 0.05
  )
  plot(0,type="n",bty="n")
  legend("topleft",legend = danio_modules_table$newname, pch = 21, pt.bg = danio_modules_table$newcolor, col = darken(danio_modules_table$newcolor, .5), ncol = 4)
  dev.off()
  return(danio_modules_table)
}

#WGCNA TF connectivity
#mks=mkss,hclust.md=hclust.method,min.md.size=30,group=libr[i]
run_WGCNA_TF= function(obj,mks,hclust.md,min.md.size,group){
  #{r setup, include=FALSE}
  
  #dir <- '/mnt/sda/alberto/projects/danio_cisreg/'
  #fcha <- function(){ gsub("-","", Sys.Date()) }
  
  knitr::opts_chunk$set(echo = TRUE)
  knitr::opts_knit$set(root.dir = dir)
  options(scipen=999)
  
  
  ## Loading Necessary Packages
  
  #{r warning = FALSE, message=FALSE}
  library(ComplexHeatmap)
  library(circlize)
  library(colorspace)
  library(ggplot2)
  library(effects)
  
  
  #Load data
  load(paste0("data/",group,"/rda/danio_wgcna_all.rda"))
  
  ## Getting cluster information for pseudobulk
  
  #scRNA-seq:
  DefaultAssay(obj)<-"RNA"
  Idents(obj)=obj$merged_sub.anno_type
  
  
  danio_tfs=data.frame(
    "id"=row.names(obj@assays$RNA),
    "class"=obj@assays$RNA@meta.features$TF  
  )
  danio_tfs=na.omit(danio_tfs)
  

  #We subset the gene expression pseudo-bulk matrix to retrieve expression from the TFs.
  
  #{r}
  danio_tfs_cw <-
    danio_counts_norm_cw[
      rownames(danio_counts_norm_cw) %in% danio_tfs$id,
    ]
  
  
  #We can browse the expression level of different TFs using this function.
  
  #{r}
  plot_tf_danio <- function(x){
    if(x %in% rownames(danio_tfs_cw)) {
      barplot(
        height=unlist(c(
          danio_tfs_cw[
            grep(
              paste("^",x,"$",sep=""),
              rownames(danio_tfs_cw),
            ),
          ]
        )),
        col = danio_ctypes$ctype_col[match(colnames(danio_tfs_cw),danio_ctypes$ctype)],
        border = "#2F2F2F",
        las=2,
        cex.names=0.7,
        main= paste(
          x,
          " (",
          danio_tfs[grep(x,danio_tfs$id),2],
          ")\n",
          sep=""
        ),
        ylab="counts per million per cluster"
      )} else {
        stop("Name not in list of TFs.")
      }
  }
  
  
  
  #As we have expression data of many transcription factors, we can visualise the global patterns of expression using heatmaps.
  #We will do so by scaling the log-transformed expression of TFs to obtain a z-score.
  
  #{r}
  danio_tfs_genecol <-
    data.frame(
      id = rownames(danio_tfs_cw),
      ctype = apply(
        danio_tfs_cw,
        1,
        highest_val # a custom function that tells which is the highest value
      )
    )
  
  danio_tfs_genecol$ctype <- 
    factor(danio_tfs_genecol$ctype,levels = colnames(danio_tfs_cw))
  danio_tfs_genecol <- 
    danio_tfs_genecol[order(danio_tfs_genecol$ctype),]
  
  danio_tfs_fc <- danio_tfs_cw[match(danio_tfs_genecol$id,rownames(danio_tfs_cw)),]
  
  
  #And using the ComplexHeatmap package:
  
  #{r, fig.height = 8, fig.width = 4}
  col_danio_tfs_expr <- colorRamp2(
    c(0:3),
    colorRampPalette(c("#f1f5ff","#b1b8c4","#2c2b46"))(4)
  )
  
  clu_ha = HeatmapAnnotation(
    name = "cell types",
    cluster = colnames(danio_tfs_fc),show_annotation_name = F,
    col = list( cluster = setNames(danio_ctypes$ctype_col[match(colnames(danio_tfs_fc),danio_ctypes$ctype)],colnames(danio_tfs_fc)))
  )
  clu_ha@anno_list$cluster@show_legend <- FALSE
  clu_ha@anno_list$cluster@label <- NULL
  
  
  
  h1 <- Heatmap(
    name="z-score",
    t(scale(t(danio_tfs_fc))),
    col = col_danio_tfs_expr,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows=FALSE,
    cluster_columns=F,
    top_annotation = clu_ha,
    #right_annotation = danio_tfs_row_anno,
    row_title=NULL
  )
  draw(h1)
  
  
  #{r}
  pdf(paste0("./figures/",group,"/WGCNA/tfs_fc_all_supp.pdf"), height = 7, width = 3)
  draw(h1)
  dev.off()
  
  
  #{r echo = FALSE}
  # do we need any of this anymore??
  # draw(danio_cor_hm)
  # draw(danio_expr_zsco_hm)
  # 
  # pdf(paste0("./figures/WGCNA/danio_TFs_heatmap.pdf"),width = 6, height = 8)
  # draw(danio_expr_zsco_hm)
  # dev.off()
  
  
  ## Analysing TFs and module connectivity
  
  #From the definition in the original WGCNA paper, the eigengene of a given module can be understood as:
  #  "The first principal component of a given module. It can be considereded a representative of the expression profiles of the genes in that given module." (slightly adapted for clarity)
  
  #For each gene, WGCNA defines a "fuzzy" measure of module membership by correlating the expression profile to that of the module eigengenes. If this value is closer to 1 it indicates that that gene is connected to many genes of that module.
  
  #We will aggregate the average expression profiles to use as eigengenes.
  
  #We can calculate the connectivity by correlating the average module expression profiles with the expression of TFs:
  
  #{r}
  
  
  tf_eigen <- 
    WGCNA::signedKME(
      scale(t(danio_tfs_cw)), # all tfs, not only those with CV > 1.25 as in wgcna markdown
      MEs, outputColumnName = ""
    )
  
  min_kme <- mks
  
  filt_top <- 
    apply(
      tf_eigen, 1,
      function(x){
        if(any(x > min_kme)){res = TRUE} else {res = FALSE} 
        return(res)
      }
    )
  
  tf_eigen <- tf_eigen[filt_top,]
  
  danio_tfs_kme <-
    data.frame(
      id = rownames(tf_eigen),
      module = apply(
        tf_eigen,
        1,
        highest_val_0 # a custom function that tells which is the highest value
      )
    )
  danio_tfs_kme$module <- factor(danio_tfs_kme$module, levels = levels(danio_id_module$module))
  
  danio_tfs_kme <- danio_tfs_kme[order(danio_tfs_kme$module),]
  
  tf_eigen <- tf_eigen[
    match(danio_tfs_kme$id,rownames(tf_eigen))
    ,
  ]
  
  #And again, we can visualise using ComplexHeatmap
  
  #{r, fig.height=8, fig.width=4, message = FALSE, warning = FALSE}
  col_kme <- 
    colorRamp2(seq(0.3,0.8,len=10),colorRampPalette(rev(viridis_pastel))(10))
  
  modules_ha <-
    HeatmapAnnotation(
      stacked = anno_barplot(
        hm_bp,
        gp = gpar(fill = danio_ctypes$ctype_col[match(colnames(hm_bp),danio_ctypes$ctype)],col=NA),
        border = FALSE,
        bar_width = 1
      ),
      show_annotation_name = F,
      annotation_name_side='right',
      gap = unit(5,"pt"),
      show_legend = FALSE
    )
  
  h2 <- Heatmap(
    name="kME",
    tf_eigen,
    col=col_kme,
    show_row_names = T,
    show_column_names = TRUE,
    cluster_rows=F,
    cluster_columns=FALSE,
    top_annotation = modules_ha,
    column_names_side = "top",
    column_names_gp = gpar(fontsize = 8)
  )
  
  draw(h2)
  
  
  
  #{r}
  pdf(paste0("./figures/",group,"/WGCNA/tfs_connectivity_all_supp.pdf"), width = 10, height = 20)
  draw(h2)
  dev.off()
  
  
  
  ## Common plot, plus TFs from the literature
  
  #We identified several TFs previously described in the literature whose region of expression within the animal is corroborated by our analyses.
  
  #{r}
  #tfs_fig2 <- read.delim2("outputs/functional_annotation/tfs_fig2.tsv", header = TRUE)
  #head(tfs_fig2)
  
  #primed gene_load
  adult_LD_DEG = vroom::vroom(paste0("./data/Table S4. List of Adult LD-DEGs.csv"))
  primed_gene_table=adult_LD_DEG %>% dplyr::filter(GC_primed == "yes")
  primed_genes=unique(primed_gene_table$zfin_id_symbol) 
  
  tfs_fig2=primed_gene_table
  
  #We will plot a conjoined heatmap of expression and connectivity, highlighting the position of these TFs in the figure.
  
  #{r, fig.height = 8, fig.width = 6, warning = FALSE, message = FALSE}
  tfs_fc_common <- t(scale(t(danio_tfs_fc)))
  tfs_fc_common <- tfs_fc_common[rownames(danio_tfs_fc) %in% rownames(tf_eigen),]
  tfs_eigen_common <- tf_eigen[rownames(tf_eigen) %in% rownames(tfs_fc_common),]
  tfs_fc_common <- tfs_fc_common[match(rownames(tfs_eigen_common),rownames(tfs_fc_common)),]
  #tfs_eigen_common <- tfs_eigen_common[match(rownames(tfs_fc_common),rownames(tfs_eigen_common)),]
  
  tfs_fig2$typecol=tfs_fig2$exp.type
  tfs_fig2$typecol[tfs_fig2$typecol=="non.spec"]="lightyellow"
  tfs_fig2$typecol[tfs_fig2$typecol=="ad.onset"]="lightblue"
  tfs_fig2$typecol[tfs_fig2$typecol=="t.tempo"]="pink"
  
  tfs_fig2=as.data.frame(tfs_fig2)
  
  
  where_tfs_common <- unlist(sapply(tfs_fig2$zfin_id_symbol,function(x){grep(x,rownames(tfs_eigen_common))}))
  tfs_rowanno_common <-
    rowAnnotation(
      TF = anno_mark(
        at = where_tfs_common,
        labels = translate_ids(names(where_tfs_common),tfs_fig2[,c(2,9)]))
    )
  
  col_danio_tfs_expr2 <- colorRamp2(
    c(0:3),
    colorRampPalette(c("white","lightyellow","firebrick"))(4)
  )
  
  
  h1_main <- Heatmap(
    name="FC",
    tfs_fc_common,
    col=col_danio_tfs_expr2,
    show_row_names = FALSE,
    show_column_names = T,
    cluster_rows=FALSE,
    cluster_columns=T,
    top_annotation = clu_ha,
    right_annotation = tfs_rowanno_common,
    bottom_annotation = clu_ha,
    column_names_side = "bottom",
    row_title=NULL
  )
  
  h2_main <- Heatmap(
    name="kME",
    tfs_eigen_common,
    col=col_kme,
    show_row_names = FALSE,
    cluster_rows=FALSE,
    cluster_columns=F,
    top_annotation = modules_ha,
    #column_labels = danio_modules_table$newname,
    column_names_side = "top",
    column_names_gp = gpar(fontsize = 7)
  )
  
  draw(h2_main+h1_main)
  # draw(h1_main+h2_main)
  
  
  
  #{r}
  pdf(paste0("./figures/",group,"/WGCNA/tfs_expr_and_connectivity_main.pdf"), 
      width = 8, height = 8)
  draw(h2_main+h1_main)
  dev.off()
  
  
  ## Saving the data
  
  #We will save the important bits for further analysis in the rest of markdowns.
  
  #{r}
  save(
    # gene expression data
    danio_tfs_cw,
    danio_tfs_fc,
    # tf data
    danio_tfs,
    # neiro_tfs,
    tfs_fig2,
    # kME
    tf_eigen,
    tfs_eigen_common,
    # visual annotations
    modules_ha,
    col_kme,
    hm_bp,
    # ctypes_rowAnno,
    # clu_ha,
    # modules_ha,
    # wg_ha,
    #pick color palette for TFs
    # destination
    file = paste0(
      "./data/",group,"/rda/tf_analysis.rda"
    )
  )
}

#WGCNA module GO
run_WGCNA_GO = function(group){
  load(paste0("data/",group,"/rda/danio_wgcna_all.rda"))
  load(paste0("./data/",group,"/rda/danio_counts.rda"))
  library(xlsx)
  
  ## Gene Ontology Analysis
  
  #For this we will use a wrapper function of the GO enrichment analysis tools provided by the package `topGO`. First the setup.
  
  #{r danio_GOs_setup, echo = FALSE, warning = FALSE}
  #gene universe
  gene_universe <- rownames(danio_counts)
  
  #danio GOdb
  library(org.Dr.eg.db)
  
  allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Dr.eg.db", ID="symbol")
  
  # gene-GO mappings
  danio_id_GO <- allGO2genes
  
  #list of genes of interest
  danio_wg_list <- split(rownames(danio_wg_module),danio_wg_module$module)
  danio_wg_list <- danio_wg_list[match(danio_modules_table$newname, names(danio_wg_list))]
  '
for (i in seq_along(danio_wg_list)) {
  file_name <- paste0( names(danio_wg_list)[i], ".txt")  # Create filename using list names
  writeLines(as.character(danio_wg_list[[i]]), paste0("./outputs/2025/wt/homer_results/RNA100_ori_modules/",file_name))  # Convert elements to character and save
}'
  
  #And now for the GO analysis:
  
  #{r danio_GOs, message = FALSE}
  # GO term analysis wrapper
  danio_wg_GO_all<- 
    getGOs(
      danio_wg_list,
      gene_universe = gene_universe,
      gene2GO = danio_id_GO
    )
  
  # gene-GO mappings
  
  #brain background
  #ref_cell: markerset1_scheir_lab (Shafer et al., 2022, Nat Ecol Evol.)
  PO_scRNA_seq=readRDS("Z:/MinK/NGS_sequencing/sc_multiome_NPO_Project_11254/old_analysis/PO_scRNA_seq_v2.rds")
  
  
  expressed_genes_HYP <- rownames(PO_scRNA_seq)[Matrix::rowSums(PO_scRNA_seq@assays$RNA$counts > 0) > 0]
  
  # GO term analysis wrapper
  danio_wg_GO_all_rosetaall_bg_elim <- 
    getGOs(
      danio_wg_list,
      gene_universe = expressed_genes_HYP,
      gene2GO = danio_id_GO,
      alg = "elim"
    )
  
  pdf(paste0("./figures/",group,"/WGCNA/wgcna_exploration_modules/danio_wgcna_module_exploration_GOs_DIFFNORM_universe_rosettaBG_elim.pdf"))
  for (i in danio_wg_GO_all_rosetaall_bg_elim$GOplot) {print(i)}
  dev.off()
  
  # GO term analysis wrapper
  danio_wg_GO_all_rosetaall_bg <- 
    getGOs(
      danio_wg_list,
      gene_universe = expressed_genes_HYP,
      gene2GO = danio_id_GO
    )
  
  pdf(paste0("./figures/",group,"/WGCNA/wgcna_exploration_modules/danio_wgcna_module_exploration_GOs_DIFFNORM_universe_rosettaBG.pdf"))
  for (i in danio_wg_GO_all_rosetaall_bg$GOplot) {print(i)}
  dev.off()
  
  # GO term analysis wrapper
  danio_wg_GO_all_danioidmodule_bg <- 
    getGOs(
      danio_wg_list,
      gene_universe = danio_id_module$id,
      gene2GO = danio_id_GO
    )
  
  pdf(paste0("./figures/",group,"/WGCNA/wgcna_exploration_modules/danio_wgcna_module_exploration_GOs_DIFFNORM_universe_wgcnaBG.pdf"))
  for (i in danio_wg_GO_all_danioidmodule_bg$GOplot) {print(i)}
  dev.off()
  
  # GO term analysis wrapper
  danio_wg_GO_all_danioidmodule_bg_elim <- 
    getGOs(
      danio_wg_list,
      gene_universe = danio_id_module$id,
      gene2GO = danio_id_GO,
      alg = "elim"
    )
  
  pdf(paste0("./figures/",group,"/WGCNA/wgcna_exploration_modules/danio_wgcna_module_exploration_GOs_DIFFNORM_universe_wgcnaBG_elim.pdf"))
  for (i in danio_wg_GO_all_danioidmodule_bg_elim$GOplot) {print(i)}
  dev.off()
  
  
  
  #Here we show a couple of GO term analysis for different modules:
  
  #{r, fig.width = 6, fig.height = 6}
  danio_wg_GO_all$GOplot$s01
  
  pdf(paste0("./figures/",group,"/WGCNA/wgcna_exploration_modules/danio_wgcna_module_exploration_GOs_DIFFNORM_universe_allgenes.pdf"))
  for (i in danio_wg_GO_all$GOplot) {print(i)}
  dev.off()
  
  
  #####
  
  
  ## Saving genome region files for motif enrichment
  
  #{r}
  #danio_promoters <- read.delim2(file = "/mnt/sda/alberto/projects/danio_cisreg/outputs/associate_peaks_genes/promoters.bed", header = FALSE)
  '
for (i in 1:length(danio_wg_list)) {
  newname <- gsub(" ", "_", names(danio_wg_list)[i])
  newname <- gsub("\\+", "POS", newname)
  newname <- gsub("\\&", "and", newname)
  tbl <- danio_promoters[danio_promoters$V4 %in% danio_wg_list[[i]],]
  write.table(
    tbl,
    file = paste0("outputs/wgcna/homer/promoters/promoters_module_",newname,".bed"),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  message("Done ",newname); rm(tbl)
}
'

## Saving everything

#{r save_GO_terms}
danio_wg_GO_table <- ldply(danio_wg_GO_all[[1]], .id="module")
write.xlsx(
  danio_wg_GO_table,
  file = paste0("./outputs/",group,"/",fcha(),"_danio_wgcna_GOterms.xlsx"),
  sheetName = "schmidtea_wgcna_GOterms",
  col.names = TRUE, row.names = FALSE, showNA = TRUE
)

pdf(paste0("./figures/",group,"/WGCNA/danio_wgcna_GOs_barplots.pdf"))
for (i in danio_wg_GO_all$GOplot) {print(i)}
dev.off()

danio_wgcna_moduleinfo <- merge(danio_id_module, expressed_genes_HYP ,by = 1,)
write.table(
  danio_wgcna_moduleinfo[danio_wgcna_moduleinfo$gene_type == "hconf",],
  file = paste0("./outputs/",group,"/danio_wgcna_modules_info.tsv"),
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)




#{r save_all_wgcna}
save(
  danio_modules_table,
  danio_id_module,
  datExpr,
  danio_wg_module,
  datKME,
  danio_wg_GO_all,
  geneTree,
  MEList,
  MEs,
  #sft,
  danio_wg_module_viz,
  danio_wg_list,
  adjacency,
  TOM,
  dissTOM,
  TOM_2,
  hm_bp,
  file = paste0("data/",group,"/rda/danio_wgcna_all.rda")
)


#{r save_subset_of_wgcna_for_graph_analysis}
save(
  danio_modules_table,
  danio_id_module,
  datExpr,
  danio_wg_module,
  datKME,
  MEs,
  danio_wg_GO_all,
  hm_bp,
  file = paste0("data/",group,"/rda/danio_wgcna.rda")
)
}

#module list export

module_out=function(group){
  #read table
  library(readxl)
  modules=read_excel(paste0("./outputs/",group,"/danio_wgcna_id_module.xlsx"),)
  mds=unique(modules$module)
  dir.create(paste0("./outputs/",group,"/wgcna_modules"))
  
  for (i in mds) {
    tb=dplyr::filter(modules, module == i)
    write.table(tb$id, file = paste0("./outputs/",group,"/wgcna_modules/",i,".txt"), sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)
  }
  #mds=paste0(vs,"_",mds)
  write.table(mds, file = paste0("./outputs/",group,"/WGCNA_module_list.txt"), sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)
}



#compare modules, hypergeometric test

#' comparemodules: gene family enrichment between modules, 
#' using the gfam object and the modules from each species, 
#' using a binomial statistical test and/or a hypergeometric test.
#' comparemodules(ma, mb, f, ga, gb) --> list(a_f, b_f, ma_f, mb_f, 
#' matrix_hypgeom, matrix_binomial, age_common, age_exclusive)
#' 
#' @param ma association file genes sp a -- gene modules
#' @param mb association file genes sp b -- gene modules
#' @param f association file for gene spp a,b -- gene family
#' 
#' ma: gene modules of species A
#' mb: gene modules of species B
#' f: Gene Family translation layer between species (GFs)
#' 
comparemodules <- function(ma,mb){
  
  
  gene_POP <- length(unique(c(ma$id,mb$id))) # population size
  
  #' Create a matrix to store stats
  fon = data.frame(
    module_a = "none", 
    module_b = "none", 
    success_in_samples = 0, 
    sample_size = 0, 
    success_in_pop = 0, 
    gene_POP = gene_POP, 
    hypgeom_pval = numeric(1), 
    hypgeom_log = numeric(1), 
    binom_pval =  numeric(1), 
    gfams_common = "none",
    gfams_excl_a = "none",
    gfams_excl_b = "none"
  )
  
  #' Transform into lists for practicality
  ma_list <- split(ma$id, ma$module) 
  mb_list <- split(mb$id, mb$module)
  
  #' Create a matrix to store results (-logpval or similar)
  PVS <- data.frame() # create matrix to store result pval
  PVS_binom <- data.frame()
  
  
  #' Compute the upper tail of the hypergeometric 
  #' distribution (survival function) for each pair of modules
  #'  as a metric of enrichment #comment from panos @ skarmetalab
  for (i in 1:length(ma_list)){ #need to dramatically optimise speed
    
    a_modulei_name <- paste0("group1",names(ma_list[i]))
    a_modulei_genes <- ma_list[[i]]
    a_fams <- unique(a_modulei_genes) #is this slow?
    
    for (j in 1:length(mb_list)){
      
      b_modulej_name <- paste0("group2",names(mb_list[j]))
      b_modulej_genes <- mb_list[[j]]
      b_fams <- unique(b_modulej_genes) #is this slow?
      
      common_fams <- paste(
        a_modulei_genes[which(a_modulei_genes %in% b_modulej_genes)], collapse = ", "
      )
      exclusive_fams_a <- paste(
        a_modulei_genes[which(!(a_modulei_genes %in% b_modulej_genes))], collapse = ", "
      )
      exclusive_fams_b <- paste(
        b_modulej_genes[which(!(b_modulej_genes %in% a_modulei_genes))], collapse = ", "
      )
      
      success_in_sample <- length(which(a_modulei_genes %in% b_modulej_genes))
      success_in_pop <- length(b_modulej_genes)
      sample_size <- length(a_modulei_genes)
      
      
      
      if (success_in_sample > 0) {
        hypg <- phyper( # HYPGEOMTEST
          # from stackoverflow/questions/8382806/hypergeometric-test-phyper
          q = success_in_sample - 1, # no. of success balls drawn from urn
          m = success_in_pop, # no. of success balls in the urn
          n = gene_POP-success_in_pop, # no. of non-success in the urn
          k = sample_size, # no. of balls drawn from urn
          lower.tail = FALSE
        )
        binom <- binom.test(
          x = success_in_sample, 
          n = sample_size, 
          p = success_in_pop / gene_POP, 
        )$p.value
        
        fon <- rbind(
          fon, 
          c(
            a_modulei_name, b_modulej_name, 
            success_in_sample, sample_size, 
            success_in_pop, gene_POP, as.numeric(hypg), 
            -log(hypg), as.numeric(binom), common_fams, # add here which are the names of the gfams enriched.
            exclusive_fams_a, # exclusive fams in a
            exclusive_fams_b # exclusive fams in b
          )
        )
      } else {
        hypg <- 1
        binom <- 1
      }
      
      PVS[i, j] <- hypg
      PVS_binom[i, j] <- binom
      
    }
  }
  
  # Tidy up of data 
  rownames(PVS) <- names(ma_list)
  colnames(PVS) <- names(mb_list)
  rownames(PVS_binom) <- names(ma_list)
  colnames(PVS_binom) <- names(mb_list)
  fon <- fon[!(fon$module_a == "none"), ]
  fon$hypgeom_pval <- as.numeric(fon$hypgeom_pval)
  fon$hypgeom_log <- as.numeric(fon$hypgeom_log)
  fon$binom_pval <- as.numeric(fon$binom_pval)
  
  loghypg <- -log(PVS)
  loghypg[loghypg > 30] <- 30
  
  logbinom <- -log(PVS_binom)
  logbinom[logbinom > 30] <- 30
  
  
  res <- list(
    stats = fon, 
    hypgeom = PVS, 
    binom = PVS_binom, 
    loghypg = loghypg, 
    logbinom = logbinom
  )
  
  return(res)
}



pars_homer=function(pth){
  
  #for homer
  #homer results folder
  homer_fd= pth
  #known motif results.
  known_HR_all=c(paste0(homer_fd,sort(list.files(homer_fd)),"/knownResults.txt"))
  
  #modules= str_split(string = sort(list.files(homer_fd)),
  #                   pattern = "_") %>% sapply(function(x) x[4])
  modules= str_split(string = sort(list.files(homer_fd)),
                     pattern = "_") %>% sapply(function(x) x[1])
  
  #load know motif files
  i=1
  tsv_motifs_all=vroom::vroom(known_HR_all[i])
  tsv_motifs_all$module = modules[i]
  tsv_motifs_all=tsv_motifs_all %>% relocate(module, .before = everything())
  colnames(tsv_motifs_all) <- 
    c(
      "module",
      "motif",
      "Consensus",
      "pvalue",
      "logpval",
      "qval",
      "no_target_seqs_with_motif",
      "pct_target_seqs_with_motif",
      "no_bg_seqs_with_motif",
      "pct_bg_seqs_with_motif"
    )
  for (i in 2:length(known_HR_all)) {
    tsv_motifs=vroom::vroom(known_HR_all[i])
    tsv_motifs$module = modules[i]
    tsv_motifs=tsv_motifs %>% relocate(module, .before = everything())
    
    colnames(tsv_motifs) <- 
      c(
        "module",
        "motif",
        "Consensus",
        "pvalue",
        "logpval",
        "qval",
        "no_target_seqs_with_motif",
        "pct_target_seqs_with_motif",
        "no_bg_seqs_with_motif",
        "pct_bg_seqs_with_motif"
      )
    tsv_motifs_all=rbind(tsv_motifs_all,tsv_motifs)
  }
  return(tsv_motifs_all)
}

###
motif_enrichdot=function(motif_res,group){

  load(file = paste0("data/",group,"/rda/danio_wgcna.rda"))
  load(file = paste0("./data/",group,"/rda/tf_analysis.rda"))
  dir.create(paste0("./figures/",group,"/WGCNA/motifs"))
  motifs_modules_prom_all=motif_res
  colnames(motifs_modules_prom_all)[1] <- "module"
  #```
  #Here is a subset from this results table, containing all the motifs with q.value values in the top 25% (==Q1) of the values distribution:
  
  #```{r}
  
if(nrow(motifs_modules_prom_all)>20){
  threshold <- quantile(motifs_modules_prom_all$qval )[3]
}else{
  threshold <- 0.1
} 
  top_motifs <- 
    motifs_modules_prom_all[
      motifs_modules_prom_all$qval  < threshold,
    ]
  
  # Transform module names column into ordered factor much like other plots from the WGCNA analysis
  top_motifs$module <- 
    factor(top_motifs$module, levels = levels(danio_modules_table$newname))
  
  # Arrange by module name and decreasing values of % peaks with motif found, to facilitate readability
  top_motifs <- 
    top_motifs[
      with(top_motifs,order(module,-as.numeric(pct_target_seqs_with_motif))),
    ]#https:/stackoverflow.com/questions/16205232/order-data-frame-by-columns-in-increasing-and-decreasing-order
  
  # Transform motifs names column into ordered factor so that ggplot respects the ordering we did in the step immediately above
  top_motifs$motif <- 
    factor(top_motifs$motif, levels = rev(unique(top_motifs$motif)))
  #```
  
  #Here the plot:
  
  #```{r, fig.width = 6, fig.height = 8}
  motifs_plot_top <- 
    ggplot(top_motifs,aes(y = motif, x = module, size = as.factor(size_pct)))+
    geom_point(color = alpha("#0b0b0b",0.5))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_size_manual(values = c(1, 2, 3, 4,5), labels = c("1-10", "10-20", "20-50", "50+"))+
    guides(size=guide_legend(title="% peaks with motif"))+
    ggtitle("Motif enrichment analysis (top motifs, qval < top 25%)")
  print(motifs_plot_top)
  #```
  
  #```{r}"
  #dir.create(paste0("./figures/motifs"))
  pdf(paste0("./figures/",group,"/WGCNA/motifs/",group,"_motifs_main.pdf"),height = 8, width = 6)
  print(motifs_plot_top)
  dev.off()
  #```
  
  #Here is the same, but for the FULL of all significantly enriched motifs (qvalue < 0.1)
  
  #```{r}
  # Transform module names column into ordered factor much like other plots from the WGCNA analysis
  motifs_modules_prom_all$module <- 
    factor(motifs_modules_prom_all$module, levels = levels(danio_modules_table$newname))
  
  # Arrange by module name and decreasing values of % peaks with motif found, to facilitate readability
  motifs_modules_prom_all <- 
    motifs_modules_prom_all[
      with(motifs_modules_prom_all,order(module,-as.numeric(pct_target_seqs_with_motif))),
    ]#https:/stackoverflow.com/questions/16205232/order-data-frame-by-columns-in-increasing-and-decreasing-order
  
  # Transform motifs names column into ordered factor so that ggplot respects the ordering we did in the step immediately above
  motifs_modules_prom_all$motif <- 
    factor(motifs_modules_prom_all$motif, levels = rev(unique(motifs_modules_prom_all$motif)))
  #```
  
  #Here is the plot
  
  #```{r, fig.height = 18, fig.width = 9}
  #motifs_modules_prom_all$target_peakset_name=factor(motifs_modules_prom_all$module,
  #                                                   levels= c("17_sst1.1_cont_pre_peak","17_avp_cont_post_peak", "17_sst1.1_cont_post_peak","17_crhb_cont_post_peak",
  #                                                             "17_avp_bPAC_pre_peak","17_sst1.1_bPAC_pre_peak","17_crhb_bPAC_pre_peak","17_avp_bPAC_post_peak"               
  #                                                                    ))
  motifs_all_plot <- 
    ggplot(
      data = motifs_modules_prom_all, 
      aes(x = module, y = motif,
          color = `pvalue`, size = as.factor(size_pct))) +
    geom_point() +
    scale_color_gradient(
      low = alpha("#401af0",0.6), high = alpha("#be0143",0.6)
    ) +
    theme_bw() +
    xlab("") + ylab("") +
    ggtitle("Motif enrichment analysis (qvalue <0.1)") +
    theme(
      text = element_text(size=10), legend.text = element_text(size=10), 
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    ) +
    scale_size_manual(values = c(1, 2, 3, 4,5), labels = c("1-10", "10-20", "20-50", "50+")) +
    guides(size=guide_legend(title="% peaks with motif"))
  
  print(motifs_all_plot)
  #```
  
  #```{r}
  # Supp Panel
  pdf(paste0("./figures/",group,"/WGCNA/motifs/",group,"motif_supp.pdf"), height = 5, width = 5)
  print(motifs_all_plot)
  dev.off()
  
  
  #####motif enrichment
  motifs <-
    motifs_modules_prom_all[
      motifs_modules_prom_all$motif %in% motifs_modules_prom_all$motif,
    ]
  #```
  
  #Here we parse and transform the table we made above table to extract two matrices, one for motif qvalues and another for percentage of peaks with enriched motif:
  
  #```{r}
  # Pivot the data to create the matrix of percentages
  matrix_pct <- as.data.frame(
    motifs[!duplicated(motifs[,c("module","motif")]),c("module","motif","pct_target_seqs_with_motif")] %>%
      pivot_wider(names_from = motif, values_from = pct_target_seqs_with_motif, values_fill = 0)
  )
  rownames(matrix_pct) <- matrix_pct$module
  matrix_pct <- matrix_pct[, -1]
  matrix_pct <- as.matrix(matrix_pct)
  
  # Pivot the data to create the matrix of qvalues
  matrix_qval <- as.data.frame(
    motifs[!duplicated(motifs[,c("module","motif")]),c("module","motif","logqval")] %>%
      pivot_wider(names_from = motif, values_from = logqval, values_fill = 0)
  )
  rownames(matrix_qval) <- matrix_qval$module
  matrix_qval <- matrix_qval[, -1]
  matrix_qval <- as.matrix(matrix_qval)
  #```
  
  #Here is how these look:
  
  #```{r}
  matrix_pct
  #```
  
  #```{r}
  matrix_qval
  #```
  
  #We will bin the values of pct into intervals to facilitate the visualisation:
  
  #```{r}
  m <- matrix(
    as.numeric( as.character(
      cut(
        x=matrix_pct,
        breaks = c(0,0.7,10,20,50,100),
        labels = c(0,1,1.5,2,2.5), 
        right = FALSE)
    )),
    ncol = ncol(matrix_pct), nrow = nrow(matrix_pct)
  )
  dimnames(m) <- dimnames(matrix_pct)
  
  m
  #```
  
  #We arrange the rows (modules) of the abridged pct matrix and the qvalue matrix to keep the same order of modules we've been doing all the time
  
  #```{r}
  m <- m[match(levels(danio_modules_table$newname),rownames(m)),]
  matrix_qval <- matrix_qval[match(levels(danio_modules_table$newname),rownames(matrix_qval)),]
  m <- m[complete.cases(m),]
  matrix_qval <- matrix_qval[complete.cases(matrix_qval),]
  #```
  
  #We also order the columns (motifs) of these matrices by highest value using a custom function. Refer to the TFs markdown, it's the same one.
  
  #```{r}
  where_highest <- data.frame(motif = rownames(t(matrix_qval)),module = apply(t(matrix_qval),1,highest_val))
  where_highest$module <- factor(where_highest$module,levels = levels(danio_id_module$module))
  where_highest <- where_highest[order(where_highest$module),]
  
  m <- m[,match(where_highest$motif,colnames(m))]
  mq <- matrix_qval[,match(where_highest$motif,colnames(matrix_qval))]
  #```
  
  #Below the different annotations to create the heatmap:
  
  #```{r}
  # Colouring function for logqvalue in heatmap
  col_fun = circlize::colorRamp2(breaks=c(2,4,6,8,10),colors=colorRampPalette(c("#401af0","#be0143"))(5))
  danio_ctypes <- read.csv(ctype_path,row.names = 1)
  
  modules_ha_motifhm <- 
    HeatmapAnnotation(
      stacked = anno_barplot(
        hm_bp[rownames(hm_bp) %in% rownames(m),],
        gp = gpar(fill = danio_ctypes$ctype_col[match(colnames(danio_wg_module),danio_ctypes$ctype)],col=NA),
        bar_width = 1
      ),
      show_legend = FALSE
    ) 
  
  # Motif box annotations (nucleotide sequences)
  motif_box_ <- motifs[motifs$motif %in% colnames(m),]
  motifs_row_anno <-
    HeatmapAnnotation(
      box = 
        anno_text(
          motif_box_$Consensus[match(colnames(m),motif_box_$motif)], # nucleotide sequences in the same order as our matrices
          gp = gpar(
            fontsize = 8, fill="#F5F5F5",
            col = translate_ids(motif_box_$module[match(colnames(m),motif_box_$motif)]
                                ,danio_modules_table[,c(4,7)]) # nucleotide sequences in the same order as our matrices
          )),
      which = "row"
    )
  #```
  
  #And the motif heatmap in the style of a dot blot:
  
  #```{r, fig.height = 8, fig.width = 8}
  motif_hm <- Heatmap(
    t(mq), # we flip it to put motifs in rows
    name = "logqvalue",
    col = col_fun,
    rect_gp = gpar(type = "none"), 
    cell_fun = function(j, i, x, y, width, height, fill) { # the embedded heatmap of dot size: basically either draw dots of certain sizes or do nothing
      if(t(m)[i,j] == 0){
        NULL
      } else{
        grid.circle(
          x = x, y = y, r = t(m)[i,j]/3 * min(unit.c(width, height)), 
          gp = gpar(fill = col_fun(t(mq)[i, j]), col = NA)
        )
      }
    },
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(cex=0.7),
    # column_labels = translate_ids(rownames(mq),danio_modules_table[,c(4,9)]),
    column_names_side = "top",
    top_annotation = modules_ha_motifhm,
    right_annotation = motifs_row_anno
  )
  draw(motif_hm,background = "transparent")
  #```
  
  #```{r}
  pdf(paste0("./figures/",group,"/WGCNA/motifs/motif_lit_main.pdf"),height = 8, width = 8)
  draw(motif_hm)
  dev.off()
  #```
  
  ### TF connectivity of literature TFs
  #primed gene_load
 ' adult_LD_DEG = vroom::vroom("./data/Table S4. List of Adult LD-DEGs.csv")
  primed_gene_table=adult_LD_DEG %>% dplyr::filter(GC_primed == "yes")
  primed_genes=unique(primed_gene_table$zfin_id_symbol) 
  load(file = paste0("./data/",group,"/rda/tf_analysis.rda"))
  tfs_fig2=primed_gene_table
  
  #We will plot a conjoined heatmap of expression and connectivity, highlighting the position of these TFs in the figure.
  
  #{r, fig.height = 8, fig.width = 6, warning = FALSE, message = FALSE}
  tfs_fc_common <- t(scale(t(danio_tfs_fc)))
  tfs_fc_common <- tfs_fc_common[rownames(danio_tfs_fc) %in% rownames(tf_eigen),]
  tfs_eigen_common <- tf_eigen[rownames(tf_eigen) %in% rownames(tfs_fc_common),]
  tfs_fc_common <- tfs_fc_common[match(rownames(tfs_eigen_common),rownames(tfs_fc_common)),]
  #tfs_eigen_common <- tfs_eigen_common[match(rownames(tfs_fc_common),rownames(tfs_eigen_common)),]
  
  tfs_fig2$typecol=tfs_fig2$exp.type
  tfs_fig2$typecol[tfs_fig2$typecol=="non.spec"]="lightyellow"
  tfs_fig2$typecol[tfs_fig2$typecol=="ad.onset"]="lightblue"
  tfs_fig2$typecol[tfs_fig2$typecol=="t.tempo"]="pink"
  
  tfs_fig2=as.data.frame(tfs_fig2)
  
  
  #```{r, fig.height = 4, fig.width = 8}
  matrix_connectivity <- tfs_eigen_common[rownames(tfs_eigen_common) %in% tfs_fig2$zfin_id_symbol,]
  
  matrix_connectivity <-
    matrix_connectivity[
      order(apply(matrix_connectivity,1,function(x){which(x==max(x))})),
    ]
  '
  ##homer_TF
  homer_TF_db= vroom::vroom("./data/homer_motifTable.txt")
  homer_TF_db$trim.Name=tidyup_motifnames(homer_TF_db$Name)
  
  #top_motif2TF
  top_motif2TF= dplyr::filter(homer_TF_db,trim.Name %in% where_highest$motif)
  
  #orthologs_info
  zf_human_ortho= vroom::vroom("./data/orthlogs_danRer11.pep__v__GRCh38.p13.pep.tsv")
  
  #top_motif2TF_danio
  zfin_id_symbol=c()
  for (i in 1:nrow(top_motif2TF)) {
    a=top_motif2TF$`Gene Symbol`[i]
    mf=unlist(str_split(a,"\\,"))
    if (length(mf) == 1) {
      x=zf_human_ortho$danRer11.pep[grepl(top_motif2TF$`Gene Symbol`[i],zf_human_ortho$GRCh38.p13.pep)]
      x=unique(
        gsub(pattern = " ",replacement = "",
        c(unlist(str_split(x,"\\||\\,")))
        )
        )
      x=paste(x,collapse = ",")
    }else{
      for (k in mf) {
        x=zf_human_ortho$danRer11.pep[grepl(k,zf_human_ortho$GRCh38.p13.pep)]
        x=gsub(pattern = " ",replacement = "",c(unlist(str_split(x,"\\||\\,"))))
        x=unique(c(x,x))
        x=paste(x,collapse = ",")
      }
    }
    zfin_id_symbol=c(zfin_id_symbol,paste(x,sep = "_"))
  }
  
  top_motif2TF$zfin_id_symbol=zfin_id_symbol
  
  #{r, fig.height = 8, fig.width = 6, warning = FALSE, message = FALSE}
  tfs_fc_common <- t(scale(t(danio_tfs_fc)))
  tfs_fc_common <- tfs_fc_common[rownames(danio_tfs_fc) %in% rownames(tf_eigen),]
  tfs_eigen_common <- tf_eigen[rownames(tf_eigen) %in% rownames(tfs_fc_common),]
  tfs_fc_common <- tfs_fc_common[match(rownames(tfs_eigen_common),rownames(tfs_fc_common)),]
  #tfs_eigen_common <- tfs_eigen_common[match(rownames(tfs_fc_common),rownames(tfs_eigen_common)),]
  
  matrix_connectivity <- tfs_eigen_common[rownames(tfs_eigen_common) %in% 
                                             unique(unlist(str_split(paste(top_motif2TF$zfin_id_symbol,collapse = ","),"\\,"))),]
  
  matrix_connectivity <-
    matrix_connectivity[
      order(apply(matrix_connectivity,1,function(x){which(x==max(x))})),
    ]
  
  
  
  #```
  
 ' 
  #```{r, fig.height = 4, fig.width = 8}
  tfs_lit_ids <- 
    setNames(
      translate_ids(rownames(matrix_connectivity),dict = tfs_fig2[,c(1,3)]),
      rownames(matrix_connectivity)
    )
  
  tf_box_rowanno <-
    HeatmapAnnotation(
      box = 
        anno_text(
          tf_box$box[match(tfs_lit_ids,tf_box$TF)],
          gp = gpar(
            fontsize = 8,
            col = tf_box$col[match(tfs_lit_ids,tf_box$TF)]
          )
        ),
      which = "row"
    )
 ' 
  kme_hm <-
    Heatmap(
      name="kME",
      matrix_connectivity,
      col=col_kme,
      show_row_names = TRUE,
      cluster_rows=T,clustering_distance_rows = "euclidean", #"euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
      cluster_columns=FALSE,
      #row_labels = tfs_lit_ids,
      top_annotation = modules_ha,
      #right_annotation = tf_box_rowanno,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 8),
      column_names_side = "top"
    )
  
  pdf(paste0("./figures/",group,"/WGCNA/motifs/kme_literature.pdf"),height = 10, width = 10)
  draw(kme_hm)
  dev.off()
  #```
  
  #And a version of the same heatmap but only with the modules which got enriched motifs (i.e. the ones we have in the motif heatmap)
  
  #```{r, fig.height = 6, fig.width = 8}
  matrix_connectivity_slim <- matrix_connectivity[,colnames(matrix_connectivity) %in% rownames(m)]
  modules_ha2_slim <- 
    HeatmapAnnotation(
      stacked = anno_barplot(
        hm_bp[rownames(hm_bp) %in% rownames(m),],
        gp = gpar(fill = danio_ctypes$ctype_col[match(colnames(danio_wg_module)[1:61],danio_ctypes$ctype)],col=NA),
        border = FALSE,
        bar_width = 1
      ),
      annotation_name_side='right',
      gap = unit(5,"pt"),
      show_legend = FALSE
    )
  
  kme_hm_slim <- Heatmap(
    name="kME",
    matrix_connectivity_slim,
    col=col_kme,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows=T,clustering_distance_rows = "euclidean", #"euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
    cluster_columns=FALSE,
    column_names_side = "top",
    column_names_rot = 0,
    #row_labels = tfs_lit_ids,
    top_annotation = modules_ha2_slim,
    #right_annotation = tf_box_rowanno,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8)
  )
  #```
  motif_hm <- Heatmap(
    t(mq), # we flip it to put motifs in rows
    name = "logqvalue",
    col = col_fun,
    rect_gp = gpar(type = "none"), 
    cell_fun = function(j, i, x, y, width, height, fill) { # the embedded heatmap of dot size: basically either draw dots of certain sizes or do nothing
      if(t(m)[i,j] == 0){
        NULL
      } else{
        grid.circle(
          x = x, y = y, r = t(m)[i,j]/3 * min(unit.c(width, height)), 
          gp = gpar(fill = col_fun(t(mq)[i, j]), col = NA)
        )
      }
    },
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(cex=0.7),
    # column_labels = translate_ids(rownames(mq),danio_modules_table[,c(4,9)]),
    column_names_side = "top",
    #top_annotation = modules_ha_motifhm,
    right_annotation = motifs_row_anno
  )  
  #```{r, fig.height = 8, fig.width = 6}
  draw(
    kme_hm_slim %v%
      motif_hm
  )
  #```
  
  #```{r, fig.height = 7, fig.width = 6}
  pdf(paste0("./figures/",group,"/WGCNA/motifs/kme_motif_literature.pdf"),height = 7, width = 6)
  draw(
    kme_hm_slim %v%
      motif_hm
  )
  dev.off()
  #```

  #```{r}
  save(
    motifs_modules_prom_all,
    motifs,
    top_motif2TF,
    matrix_pct,
    matrix_qval,
    file = paste0("data/",group,"/rda/motif_analysis.rda")
  )

}



TF_primed= function(group){
  #{r warning = FALSE, message=FALSE}
  library(ComplexHeatmap)
  library(circlize)
  library(colorspace)
  library(ggplot2)
  library(effects)
  
  
  #Load data
  #TFs
  danio_tfs=data.frame(
    "id"=row.names(ss1@assays$RNA),
    "class"=ss1@assays$RNA@meta.features$TF  
  )
  danio_tfs=unique(na.omit(danio_tfs))
  
  load(file = paste0("data/",group,"/rda/danio_wgcna.rda"))
  load(paste0("data/",group,"/rda/danio_wgcna_all.rda"))
  load(paste0("./data/",group,"/rda/danio_counts.rda"))
  
  #We subset the gene expression pseudo-bulk matrix to retrieve expression from the TFs.
  
  #{r}
  danio_tfs_cw <-
    danio_counts_norm_cw[
      rownames(danio_counts_norm_cw) %in% danio_tfs$id,
    ]
  
  
  #We can browse the expression level of different TFs using this function.
  
  #{r}
  plot_tf_danio <- function(x){
    if(x %in% rownames(danio_tfs_cw)) {
      barplot(
        height=unlist(c(
          danio_tfs_cw[
            grep(
              paste("^",x,"$",sep=""),
              rownames(danio_tfs_cw),
            ),
          ]
        )),
        col = danio_ctypes$simp_ctype_col[match(colnames(danio_tfs_cw),danio_ctypes$np_sub_ctype)],
        border = "#2F2F2F",
        las=2,
        cex.names=0.7,
        main= paste(
          x,
          " (",
          danio_tfs[grep(x,danio_tfs$id),2],
          ")\n",
          sep=""
        ),
        ylab="counts per million per cluster"
      )} else {
        stop("Name not in list of TFs.")
      }
  }
  
  
  
  #As we have expression data of many transcription factors, we can visualise the global patterns of expression using heatmaps.
  #We will do so by scaling the log-transformed expression of TFs to obtain a z-score.
  
  #{r}
  danio_tfs_genecol <-
    data.frame(
      id = rownames(danio_tfs_cw),
      ctype = apply(
        danio_tfs_cw,
        1,
        highest_val # a custom function that tells which is the highest value
      )
    )
  
  danio_tfs_genecol$ctype <- 
    factor(danio_tfs_genecol$ctype,levels = colnames(danio_tfs_cw))
  danio_tfs_genecol <- 
    danio_tfs_genecol[order(danio_tfs_genecol$ctype),]
  
  danio_tfs_fc <- danio_tfs_cw[match(danio_tfs_genecol$id,rownames(danio_tfs_cw)),]
  
  
  #And using the ComplexHeatmap package:
  
  #{r, fig.height = 8, fig.width = 4}
  col_danio_tfs_expr <- colorRamp2(
    c(0:3),
    colorRampPalette(c("#f1f5ff","#b1b8c4","#2c2b46"))(4)
  )
  
  clu_ha = HeatmapAnnotation(
    name = "cell types",
    cluster = colnames(danio_tfs_fc),show_annotation_name = F,
    col = list( cluster = setNames(danio_ctypes$ctype_col[match(colnames(danio_tfs_fc),danio_ctypes$ctype)],colnames(danio_tfs_fc)))
  )
  clu_ha@anno_list$cluster@show_legend <- FALSE
  clu_ha@anno_list$cluster@label <- NULL
  
  
  
  h1 <- Heatmap(
    name="z-score",
    t(scale(t(danio_tfs_fc))),
    col = col_danio_tfs_expr,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows=FALSE,
    cluster_columns=F,
    top_annotation = clu_ha,
    #right_annotation = danio_tfs_row_anno,
    row_title=NULL
  )
  draw(h1)
  
  
  #{r}
  pdf(paste0("./figures/",group,"/WGCNA/tfs_fc_all_supp.pdf"), height = 7, width = 3)
  draw(h1)
  dev.off()
  
  
  #{r echo = FALSE}
  # do we need any of this anymore??
  # draw(danio_cor_hm)
  # draw(danio_expr_zsco_hm)
  # 
  # pdf(paste0("./figures/WGCNA/danio_TFs_heatmap.pdf"),width = 6, height = 8)
  # draw(danio_expr_zsco_hm)
  # dev.off()
  
  
  ## Analysing TFs and module connectivity
  
  #From the definition in the original WGCNA paper, the eigengene of a given module can be understood as:
  #  "The first principal component of a given module. It can be considereded a representative of the expression profiles of the genes in that given module." (slightly adapted for clarity)
  
  #For each gene, WGCNA defines a "fuzzy" measure of module membership by correlating the expression profile to that of the module eigengenes. If this value is closer to 1 it indicates that that gene is connected to many genes of that module.
  
  #We will aggregate the average expression profiles to use as eigengenes.
  
  #We can calculate the connectivity by correlating the average module expression profiles with the expression of TFs:
  
  #{r}
  
  tf_eigen <- 
    WGCNA::signedKME(
      scale(t(danio_tfs_cw)), # all tfs, not only those with CV > 1.25 as in wgcna markdown
      MEs, outputColumnName = ""
    )
  
  min_kme <- .55
  
  filt_top <- 
    apply(
      tf_eigen, 1,
      function(x){
        if(any(x > min_kme)){res = TRUE} else {res = FALSE} 
        return(res)
      }
    )
  
  tf_eigen <- tf_eigen[filt_top,]
  
  danio_tfs_kme <-
    data.frame(
      id = rownames(tf_eigen),
      module = apply(
        tf_eigen,
        1,
        highest_val_0 # a custom function that tells which is the highest value
      )
    )
  danio_tfs_kme$module <- factor(danio_tfs_kme$module, levels = levels(danio_id_module$module))
  
  danio_tfs_kme <- danio_tfs_kme[order(danio_tfs_kme$module),]
  
  tf_eigen <- tf_eigen[
    match(danio_tfs_kme$id,rownames(tf_eigen))
    ,
  ]
  
  #And again, we can visualise using ComplexHeatmap
  
  #{r, fig.height=8, fig.width=4, message = FALSE, warning = FALSE}
  col_kme <- 
    colorRamp2(seq(0.3,0.8,len=10),colorRampPalette(rev(viridis_pastel))(10))
  
  modules_ha <-
    HeatmapAnnotation(
      stacked = anno_barplot(
        hm_bp,
        gp = gpar(fill = danio_ctypes$simp_ctype_col[match(colnames(hm_bp),danio_ctypes$ctype)],col=NA),
        border = FALSE,
        bar_width = 1
      ),
      show_annotation_name = F,
      annotation_name_side='right',
      gap = unit(5,"pt"),
      show_legend = FALSE
    )
  
  h2 <- Heatmap(
    name="kME",
    tf_eigen,
    col=col_kme,
    show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    top_annotation = modules_ha,
    column_names_side = "top",
    column_names_gp = gpar(fontsize = 8)
  )
  
  draw(h2)
  
  
  
  #{r}
  pdf(paste0("./figures/",group,"/WGCNA/tfs_connectivity_all_supp.pdf"), width = 10, height = 10)
  draw(h2)
  dev.off()
  
  
  
  ## Common plot, plus TFs from the literature
  
  #We identified several TFs previously described in the literature whose region of expression within the animal is corroborated by our analyses.
  
  #{r}
  #tfs_fig2 <- read.delim2("outputs/functional_annotation/tfs_fig2.tsv", header = TRUE)
  #head(tfs_fig2)
  
  #primed gene_load
  adult_LD_DEG = vroom::vroom("./data/Table S4. List of Adult LD-DEGs.csv")
  primed_gene_table=adult_LD_DEG %>% dplyr::filter(GC_primed == "yes")
  primed_genes=unique(primed_gene_table$zfin_id_symbol) 
  
  tfs_fig2=primed_gene_table
  
  #We will plot a conjoined heatmap of expression and connectivity, highlighting the position of these TFs in the figure.
  
  #{r, fig.height = 8, fig.width = 6, warning = FALSE, message = FALSE}
  tfs_fc_common <- t(scale(t(danio_tfs_fc)))
  tfs_fc_common <- tfs_fc_common[rownames(danio_tfs_fc) %in% rownames(tf_eigen),]
  tfs_eigen_common <- tf_eigen[rownames(tf_eigen) %in% rownames(tfs_fc_common),]
  tfs_fc_common <- tfs_fc_common[match(rownames(tfs_eigen_common),rownames(tfs_fc_common)),]
  #tfs_eigen_common <- tfs_eigen_common[match(rownames(tfs_fc_common),rownames(tfs_eigen_common)),]
  
  tfs_fig2$typecol=tfs_fig2$exp.type
  tfs_fig2$typecol[tfs_fig2$typecol=="non.spec"]="lightyellow"
  tfs_fig2$typecol[tfs_fig2$typecol=="ad.onset"]="lightblue"
  tfs_fig2$typecol[tfs_fig2$typecol=="t.tempo"]="pink"
  
  tfs_fig2=as.data.frame(tfs_fig2)
  
  
  where_tfs_common <- unlist(sapply(tfs_fig2$zfin_id_symbol,function(x){grep(x,rownames(tfs_eigen_common))}))
  tfs_rowanno_common <-
    rowAnnotation(
      TF = anno_mark(
        at = where_tfs_common,
        labels = translate_ids(names(where_tfs_common),tfs_fig2[,c(2,9)]))
    )
  
  col_danio_tfs_expr2 <- colorRamp2(
    c(0:3),
    colorRampPalette(c("white","lightyellow","firebrick"))(4)
  )
  
  
  h1_main <- Heatmap(
    name="FC",
    tfs_fc_common,
    col=col_danio_tfs_expr2,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows=FALSE,
    cluster_columns=T,
    top_annotation = clu_ha,
    right_annotation = tfs_rowanno_common,
    bottom_annotation = clu_ha,
    column_names_side = "bottom",
    row_title=NULL
  )
  
  h2_main <- Heatmap(
    name="kME",
    tfs_eigen_common,
    col=col_kme,
    show_row_names = FALSE,
    cluster_rows=FALSE,
    cluster_columns=F,
    top_annotation = modules_ha,
    # column_labels = danio_modules_table$codename,
    column_names_side = "top",
    column_names_gp = gpar(fontsize = 7)
  )
  
  draw(h2_main+h1_main)
  # draw(h1_main+h2_main)
  
  
  
  #{r}
  pdf(paste0("./figures/",group,"/WGCNA/tfs_expr_and_connectivity_main.pdf"), width = 10, height = 8)
  draw(h2_main+h1_main)
  dev.off()
  
  
  ## Saving the data
  
  #We will save the important bits for further analysis in the rest of markdowns.
  
  #{r}
  save(
    # gene expression data
    danio_tfs_cw,
    danio_tfs_fc,
    # tf data
    danio_tfs,
    # neiro_tfs,
    tfs_fig2,
    # kME
    tf_eigen,
    tfs_eigen_common,
    # visual annotations
    modules_ha,
    col_kme,
    hm_bp,
    # ctypes_rowAnno,
    # clu_ha,
    # modules_ha,
    # wg_ha,
    #pick color palette for TFs
    # destination
    file = paste0(
      "./data/",group,"/rda/",
      "tf_analysis.rda"
    )
  )
}


