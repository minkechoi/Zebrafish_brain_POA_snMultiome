# =============================================================================
# 08_a_wgcna_prep_psudobilk.R
# -----------------------------------------------------------------------------
# Purpose : Prepare per-library pseudobulk input for WGCNA (step 08). Defines
#           wgcna_prep(): for one library it pseudobulks RNA (and ATAC) counts by
#           cell type, computes cell-specificity weights, DESeq2-normalises,
#           assigns each gene to its highest-expressing cell type, builds a
#           cell-type co-occurrence tree/heatmap, and saves everything for WGCNA.
# Function: wgcna_prep(obj, group, ctype_info_path, mxtype, hclust.md)
#   obj             - Seurat object (all libraries)
#   group           - one orig.ident to subset (a single library)
#   ctype_info_path - CSV of cell-type info (colours/order)
#   mxtype          - name of the count matrix variable to use for gene->type
#   hclust.md       - hclust linkage method (e.g. "ward.D2")
# Output  : figures/tables under ./figures/<group>/, ./data/<group>/, and the
#           bundle ./data/<group>/rda/danio_counts.rda used by 08_WGCNA.
# =============================================================================

wgcna_prep=function(obj,group,ctype_info_path,mxtype,hclust.md){

  # Subset to a single library and create its output folder scaffold
  obj= subset(obj, orig.ident == group)
  dir.create(paste0("./figures/",group,"/"))
  dir.create(paste0("./figures/",group,"/WGCNA/"))
  dir.create(paste0("./figures/",group,"/WGCNA/wgcna_exploration_modules/"))
  dir.create(paste0("./figures/",group,"/WGCNA/wgcna_exploration_modules/featplots_modules_deseq_CW/"))
  dir.create(paste0("./outputs/",group,"/"))
  dir.create(paste0("./data/",group,"/"))
  dir.create(paste0("./data/",group,"/rda/"))
  dir.create(paste0("./data/",group,"/rds/"))

  ## Getting cluster information for pseudobulk

  #scRNA-seq: sum RNA counts per annotated cell type
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
  # Cells expressing each gene per cluster (for weighting) + cluster sizes
  ss_psbulk_ncells <-
    pseudobulk_ncells(
      x = obj@assays$RNA@counts,
      identities = Idents(obj),
      min_counts = 1
    )

  cluster_size =
    c(table(obj$merged_sub.anno_type))


  #scATAC-seq: same pseudobulk for the ATAC assay
  ss_atac_clusters <-
    setNames(
      obj$merged_sub.anno_type,
      colnames(obj)
    )

  ss_atac_pseudobulk <-
    pseudobulk(
      x = obj@assays$ATAC$counts,
      ident = ss_atac_clusters
    )
  ss_atac_pseudobulk=ss_atac_pseudobulk[,str_sort(colnames(ss_atac_pseudobulk), numeric = TRUE)]

  #reorder column: map each cell type to its display colour
  coltb=data.frame("id"=obj$merged_sub.anno_type,
                   "col"=obj$ctype_col)
  coltb=arrange(coltb,id)
  coltb=unique(coltb)
  #barclaobj=sapply(str_split(colnames(ss_pseudobulk),"_"), `[`, 1)

  cols=translate_ids(coltb$id,dict = coltb[,c(1,2)])

  # Barplot: number of genes quantified (>5 counts) per cell type
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

  ## Save Data (pseudobulk matrices as TSVs)


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


  write.table(
    ss_atac_pseudobulk,
    file = paste0("./data/",group,"/ss_atac_pseudobulk.tsv"),
    sep = "\t", dec = ".",
    row.names = TRUE, quote = FALSE
  )

  # 09 ----------------------------------------------------------------------
  # Normalisation + weighting block (libraries loaded here for standalone use)
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


  ## Loading cell type information (colours/order aligned to this object)


  danio_ctypes <- read.csv(ctype_info_path,row.names = 1)
  danio_ctypes$ctype <- factor(danio_ctypes$ctype,levels = unique(danio_ctypes$ctype))
  danio_ctypes=danio_ctypes[c(rownames(obj[["merged_sub.anno_type"]])),]
  danio_ctypes$merged_sub.anno_type=obj$merged_sub.anno_type
  danio_ctypes$ctype <- factor(danio_ctypes$merged_sub.anno_type,levels = unique(danio_ctypes$merged_sub.anno_type))

  ## Loading counts (drop near-zero genes)


  danio_counts <- ss_pseudobulk
  danio_counts <- danio_counts[rowSums(danio_counts) > 10,]


  #And if we check the distribution of counts per cluster:


  cisreg_psbulk_ncells <- ss_psbulk_ncells


  #We will calculate a cell weight matrix to weigh the expression values from the pseudobulk count matrix. For this it will calculate how many cells (in percentage) are expressing a given gene in a given cluster, in relation to how many cells (in percentage) are expressing that gene in the rest of clusters. We set the minimum number of total counts as 30 and the minimum cells to take into account as 1.


  min_counts <- 30
  min_cells <- 3
  cluster_size = cluster_size


  # Gene-by-cluster specificity weights
  #Now we create the weight values matrix:


  danio_psbulk_cellweights <-
    get_cellweight_matrix(
      x = as.matrix(danio_counts),
      y = cisreg_psbulk_ncells,
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
  # Boxplots of count normalisation (raw / DESeq2 / weighted / weighted+DESeq2)

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
  danio_counts_norm_cw=get(mxtype)     # choose which matrix drives the assignment

  danio_genecolor <- data.frame(
    id = rownames(danio_counts_norm_cw),
    genecolor = apply(
      danio_counts_norm_cw,
      1,
      function(x) {
        a <- which( x == max(x) )       # cell type with highest expression
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

  # set fixed seed (deterministic clustering; no bootstrap here)
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
  danio_cisreg_cooc_hm <- Heatmap(
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

  draw(danio_cisreg_cooc_hm)

  tiff(paste0("./figures/",group,"/WGCNA/cooccurrence.tiff"),
       width = 45,height = 40,units = "cm", res = 300,compression = "lzw")
  draw(danio_cisreg_cooc_hm)

  dev.off()

  ## Save the data (everything WGCNA step 08 needs for this library)


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
