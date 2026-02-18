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

DefaultAssay(ss)="RNA"

Idents(ss)=ss$merged_sub.anno_type
ss_clusters <- 
  setNames(
    ss$merged_sub.anno_type,
    colnames(ss)
  )

ss_pseudobulk <- 
  pseudobulk(
    x = ss@assays$RNA$counts,
    ident = ss_clusters
  )

ss_psbulk_ncells <- 
  pseudobulk_ncells(
    x = ss@assays$RNA@counts,
    identities = Idents(ss),
    min_counts = 1
  )


cluster_size = 
  c(table(ss$merged_sub.anno_type))

###
type_table_m.v2= read.csv(paste0("./outputs/cell_type/cell_type_table_m_modified.csv"), header = T, row.names = 1)


ss_pseudobulk

ss_pseudobulk <- ss_pseudobulk[rowSums(ss_pseudobulk) > 10,]

#And if we check the distribution of counts per cluster:
ss_psbulk_ncells 

#We will calculate a cell weight matrix to weigh the expression values from the pseudobulk count matrix. For this it will calculate how many cells (in percentage) are expressing a given gene in a given cluster, in relation to how many cells (in percentage) are expressing that gene in the rest of clusters. We set the minimum number of total counts as 30 and the minimum cells to take into account as 1.


min_counts <- 30
min_cells <- 5
cluster_size 


#Now we create the weight values matrix:
ss_psbulk_cellweights <- 
  get_cellweight_matrix(
    x = as.matrix(ss_pseudobulk),
    y = ss_psbulk_ncells,
    C = cluster_size,
    min_counts = min_counts,
    min_cells = min_cells
  )


#A quick look at this matrix:


ss_psbulk_cellweights[1:5,1:5]

dim(ss_psbulk_cellweights)


#We subset the counts matrix to keep the same genes that were retrieved in the cell weight values matrix
colnames(ss_pseudobulk)=gsub("\\/",".",colnames(ss_pseudobulk))
type_table_m.v2$numbered=gsub("\\/",".",type_table_m.v2$numbered)


m_ <- as.matrix(ss_pseudobulk)
m_ <- m_[rowSums(m_) >= min_counts,]
m_ <- m_[rownames(m_) %in% rownames(ss_psbulk_cellweights),]


#And we normalise using DESeq2


m_dds <- DESeqDataSetFromMatrix(
  countData = m_,
  colData = data.frame(condition = colnames(m_)),
  design = ~ condition)
m_dds <- estimateSizeFactors(m_dds)
ss_pseudobulk_norm <- counts(m_dds, normalized=TRUE)


#The final matrix is the log-transformed of these DESeq2-normalised values, multiplied by the cell weights.
ss_pseudobulk_norm_cw <- (log1p(ss_pseudobulk_norm) * ss_psbulk_cellweights)
#ss_pseudobulk_norm_cw <- log1p(ss_pseudobulk_norm)

# comparABle function tidyup from source
source("./ext_code/comparABle-main/comparABle-main/code/functions/tidyup_functions.R")

danio_cpm_cooc <-
  tidyup(
    ss_pseudobulk_norm_cw,#[rownames(danio_cpm) %in% danio_hvgs,],
    highlyvariable = TRUE #FALSE # TRUE if not using the subset of danio_hvgs
  )

# set fixed seed
set.seed(4343)
h <- c(0.75,0.9)
clustering_algorithm <- "hclust" #"nj", "hclust"
clustering_method <-"ward.D2"   # "complete", "average", "ward.D2"
cor_method <- "pearson"
p <- 0.01
danio_cpm_vargenes = rownames(danio_cpm_cooc)

# Levy et al 2021 'treeFromEnsembleClustering' from source
source("./ext_code/r_code/functions/treeFromEnsembleClustering.R")
library(ggtree)
set.seed(4343)
cooc <- treeFromEnsembleClustering(
  x=danio_cpm_cooc, p=p, h=h, n = 10000, vargenes = danio_cpm_vargenes, bootstrap=T,
  clustering_algorithm=clustering_algorithm, clustering_method=clustering_method, 
  cor_method=cor_method
)


tree=cooc$tree
tree$tip.label=type_table_m.v2$numbered_mk
#tree$edge.length=log1p(tree$edge.length)

node_col=as.list(type_table_m.v2$clcol)
names(node_col)=type_table_m.v2$numbered_mk
#node_col=node_col[- grep("UnD",type_table_m.v2$numbered)]

und=c(list("und"=tree$tip.label[grep("UnD",tree$tip.label)]),
         as.list(tree$tip.label[- grep("UnD",tree$tip.label)],))
names(und)=c("und",tree$tip.label[- grep("UnD",tree$tip.label)])

#ggtree(tree, branch.length="none")+
ggtree(tree)+
  geom_tiplab(align=TRUE, linetype='dashed', linesize=.3,hjust = -0.1)+
  geom_text(aes(label=node), hjust=1)

tree=groupOTU(tree,und,"group")

a= ggtree(tree) %>% ggtree::rotate(58) %>% ggtree::rotate(59)+ #flip(90,89)+ 
  geom_tiplab(align=TRUE, linetype='dashed', linesize=.3,hjust = -0.1)+
  #geom_hilight(node=c(65,68), fill="gold",colour="lightgray", alpha =0.2, to.bottom = T) + #neuron
  #geom_hilight(node=c(55,56), fill="pink",colour="lightgray", to.bottom = T)+ #olig
  #geom_hilight(node=9, fill="blue",colour="lightgray", to.bottom = T)+ #microglia
  #geom_hilight(node=125, fill="forestgreen",colour="lightgray", alpha =0.2, to.bottom = T)+ #RG/Ependymal
  #geom_hilight(node=1, fill="purple",colour="lightgray", alpha =0.2, to.bottom = T)+ #astro
  geom_tippoint(aes(color=group), size=3)+ 
  scale_color_manual(values = c(node_col)) +
  theme_tree(legend.position='none')
a


pdf(paste0("./figures/cell_type/psudobulk_hcluster_complete.pdf"), width = 1.5, height = 10)

#par(mar = c(bottom, left, top, right)) 
#par(mar = c(0.5, 0.5, 0.5, 10)) 
plot(a)
dev.off()
#with label
s= a+
  geom_cladelabel(node=1, label=".", 
                  color="white", offset=25, offset.text = 0.2,align=TRUE)
s  

taxa_order <- get_taxa_name(s)

tiff(paste0("./figures/cell_type/psudobulk_hcluster_complete_label.tiff"),
     width = 28,height = 30,units = "cm", res = 300,compression = "lzw",bg="transparent")
#par(mar = c(bottom, left, top, right)) 
#par(mar = c(0.5, 0.5, 0.5, 10)) 
plot(s)

dev.off()

pdf(paste0("./figures/cell_type/psudobulk_hcluster_complete_label.pdf"), width = 20, height = 10)
plot(s)

dev.off()

#The resulting heatmap of similarity:

#{r fig.width=7.5, fig.height=6.5, echo = FALSE}

clu_ha = HeatmapAnnotation(
  name = "celltypes",
  cluster = factor(colnames(cooc$cooccurrence), levels = unique(type_table_m.v2$numbered)),
  col = list( cluster = setNames(c(unlist(node_col),"grey"),c(translate_ids(unlist(node_col),type_table_m.v2[,c(12,8)]),"51_UnD"))),
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
  #left_annotation = clu_ha,
  top_annotation = clu_ha,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
)

draw(danio_cisreg_cooc_hm)

decorate_column_dend("co-occurence", {
  vp = current.viewport()
  yscale = vp$yscale
  grid.yaxis(at = yscale[2] - 0:10, label = 0:10)
})

tiff(paste0("./figures/cell_type/psudobulk_cooccurrence_complete.tiff"),
     width = 45,height = 40,units = "cm", res = 300,compression = "lzw")
draw(danio_cisreg_cooc_hm)

dev.off()


save(
  cooc,
  taxa_order,
  file=paste0("./data/rda/hcluster.rda")
)