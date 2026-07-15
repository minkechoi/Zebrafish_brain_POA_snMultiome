
# =============================================================================
# ananse_graph for average_v2.R
# -----------------------------------------------------------------------------
# Purpose : Same ANANSE gene-regulatory-network analysis as ananse_graph_v3.R,
#           but built on the AVERAGED network across conditions/replicates rather
#           than per-cell-type networks. Imports the ANANSE tables, overlays DEGs,
#           finds DEG-driving TFs (direct + secondary), and renders TF->target
#           heatmaps and influence figures. Uses top_edges_per_tgtf() from
#           ananse_graph_function.R.
# Inputs  : ANANSE outputs under ./scANANSE/... and DEG results (step 07/09).
# Outputs : averaged-GRN figures and TF/target tables under ./figures/, ./outputs/.
# Note    : v2 (averaged) counterpart of ananse_graph_v3.R (per-cell-type).
# =============================================================================

library(circlize)
library(colorspace)
library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(ggplot2)
library(harmony)
library(igraph)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(Seurat)
library(tidyverse)
library(topGO)
library(viridis)
library(WGCNA)
library(xlsx)
library(colorspace)

## Load libraries

#set the data and annotations 
set.seed(1234)
options(future.globals.maxSize = 120000 * 1024^2)

r.variable=4000
vs="06_2025_v4"
umap="wnn.umap"

#set working dir
if(getwd()!=paste0("D:/projects/scMultiome_oxt/",vs)){
  setwd(paste0("./",vs))
}

#functions adopted and modified from Alberto Perez-Posada @apposada
source("./ext_code/r_code/functions/sourcefolder.R")
sourceFolder("./ext_code/r_code/functions/")

#other temp function
'%!in%' <- function(x,y)!('%in%'(x,y))

ss1=readRDS(paste0("./data/rds/Step7_var",r.variable,".rds"))

###
danio_tfs <- read.delim2("./scANANSE/analysis_11_2025/gimme/danRer11.gimme.vertebrate.v5.0.motif2factors.txt")

danio_all_tf_ananse <- unique(danio_tfs$Factor)
tfs_all <- rownames(ss1@assays$RNA)[which(ss1@assays$RNA@meta.features$TF != is.na(ss1@assays$RNA@meta.features$TF))]
tfs_all=data.frame("gene"=tfs_all,"Symbol"=tfs_all)


##load()
#dir.create(paste0("./outputs/ananse/"))
#dir.create(paste0("./outputs/ananse/network/"))
nw_dir <- paste0("./scANANSE/analysis_11_2025/network/")
nw_paths <- 
  paste0(
    nw_dir,
    list.files(path = nw_dir, pattern = ".tsv")
  )
nw_paths=nw_paths[grepl("15|35.|45|average",nw_paths)]
nw_names <- gsub(".*/","",gsub(".tsv","",nw_paths))
nw_names <- nw_names[grepl("15|35.|45|average",nw_names)]

lg <- lapply(1:length(nw_names), function(x) list())
names(lg) <- nw_names


###net import : read the averaged ANANSE TF->target network (edges + influence weights)
prob_thresh <- .8

for (i in 1:length(lg)){
  message("Loading netwowk ", i)
  nw <- data.table::fread(nw_paths[i], header = TRUE)
  nw <- as.data.frame(nw)
  
  message("Subsetting network ", i)
  nw <- nw[nw$prob > prob_thresh,]
  nw$tf <- sub("—.*","",nw$tf_target)
  nw$tg <- sub(".*—","",nw$tf_target)
  
  expr <- unique(rbind(
    unique(data.frame(
      gene = nw$tf,
      expr = nw$tf_expression
    )),
    unique(data.frame(
      gene = nw$tg,
      expr = nw$target_expression
    ))
  ))
  # break ties of having ranked tfs independently in TF and TG columns, using mean
  are_dup <- which(duplicated(expr$gene))
  dups <- expr$gene[are_dup]
  d2 <- stats::aggregate(expr[expr$gene %in% dups,2],by = list(gene = expr[expr$gene %in% dups,1]), FUN = mean)
  colnames(d2) <- c("gene","expr")
  expr <- rbind(
    expr[!(expr$gene %in% dups),],
    d2
  )
  
  # we stored expression elsewhere, we remove it now
  nw <- nw[,c(7,8,2,5,6)]
  
  message("Creating graph ", i)
  g <- graph.data.frame(d = nw, directed = TRUE)
  g <- igraph::delete.vertices(g, which(igraph::degree(g) == 0))
  
  message("adding expression, graph ", i)
  V(g)$expression <- expr$expr[match(V(g)$name, expr$gene)]
  
  message("Calculating centralities and degrees of ", i, " graph")
  V(g)$centr <- relativise(nan_to_zero(closeness(g, mode = "all")))
  V(g)$outcentr_score <- relativise(nan_to_zero(closeness(g, mode = "out")))
  
  V(g)$outdegree <- nan_to_zero(igraph::degree(g, mode = "out"))
  V(g)$indegree <- nan_to_zero(igraph::degree(g, mode = "in"))
  V(g)$rel_outdegree <- V(g)$outdegree/(V(g)$outdegree+V(g)$indegree)
  
  lg[[i]] <- g
  
  rm(nw,g,expr,d2,are_dup,dups)
}


str(lg,max.level = 1)

#####color set

ref_colors <- c("#368236","#0e6655", "#229954", "#138d75")

lighten_steps <- function(color, n = 4) {
  lighten_seq <- seq(0, 0.6, length.out = n)  
  lighten(color, amount = lighten_seq)
}

color_list <- lapply(ref_colors, lighten_steps)

all_colors <- unlist(color_list)


broadcols <- c(all_colors,"grey")


names(broadcols) <- c(nw_names)

####

pdf(paste0("./figures/ananse_graphs_no_tfs_average.pdf"), height = 3.5, width = 5)
barplot(
  sapply(
    lg,
    function(x){length(V(x)$rel_outdegree[V(x)$rel_outdegree>0])}), # no. of TFs
  col = broadcols,
  border = darken(broadcols, .4),
  ylim = c(0,1000),
  ylab = "no. of TFs in graph",
  las = 2
)
dev.off()
pdf(paste0("./figures/ananse_graphs_no_tfs23_average.pdf"), height = 3.5, width = 6)

barplot(
  sapply(
    lg,
    function(x){
      length(V(x)$rel_outdegree[V(x)$rel_outdegree > .9])/
        vcount(x)
    }),
  col = broadcols,
  border = darken(broadcols, .4)
)
dev.off()


####
tf_cen <- 
  lapply(
    lg,
    function(x){
      setNames(V(x)$centr[V(x)$outdegree>0], V(x)$name[V(x)$outdegree>0])
    }
  )

tf_cen <- t(ldply(tf_cen,function(s){t(data.frame(unlist(s)))}))
colnames(tf_cen) <- tf_cen[1,]
tf_cen <- as.data.frame(tf_cen[-1,])
tf_cen = data.frame(lapply(tf_cen, function(x) as.numeric(x)),
                    check.names=F, row.names = rownames(tf_cen))
tf_cen[is.na(tf_cen)] = 0

#mt_ori=danio_tfs_pre[which(danio_tfs_pre$Factor %in% rownames(tf_cen)) ,1]
#danio_tfs
tf_cen <- tf_cen[rownames(tf_cen) %in% tfs_all$gene,]
tf_cen_rel <- t(apply(tf_cen,1,relativise))
tf_cen_rel <-
  tf_cen_rel[
    match(
      names(sort(sapply(as.data.frame(t(tf_cen_rel)),function(x){which(x == max(x))[1]}))), #staircase sorting very quick and dirty
      rownames(tf_cen_rel)
    ),
  ]

#
tf_cen_hm <- 
  Heatmap(
    name = "scaled\ngraph\ncentrality",
    tf_cen_rel,
    col = colorRampPalette(c("#f1f5ff","#b1b8c4","#2c2b46"))(10), # sequential_hcl(10, "Sunset"),
    show_row_names = T,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_names_side = "top",
    show_column_names = T#,
    #top_annotation = HeatmapAnnotation(cluster = colnames(tf_cen_rel), col = list(cluster = setNames(broadcols,colnames(tf_cen_rel))))
  )

pdf(paste0("./figures/ananse_cen_tf_hm2_average.pdf"), width = 8, height = 16)
draw(tf_cen_hm)
dev.off()
###

clu_method = "ward.D2"
tf_cen_cor_hm <- 
  Heatmap(
    name = "Pearson\n(TF centrality)",
    cor(tf_cen),
    col = colorRamp2(seq(0.4,1,length = 10),sequential_hcl(10, "Sunset")),
    clustering_method_rows = clu_method,
    clustering_method_columns = clu_method ,
    top_annotation = 
      HeatmapAnnotation(
        cluster = colnames(tf_cen_rel),
        col = list(cluster = setNames(broadcols,colnames(tf_cen_rel))),
        show_legend = FALSE
      ),
    left_annotation = 
      HeatmapAnnotation(
        cluster = colnames(tf_cen_rel),
        col = list(cluster = setNames(broadcols,colnames(tf_cen_rel))),
        which = "row"
      )
  )

pdf(paste0("./figures/ananse_cen_cor_hm2_average.pdf"), width = 7, height = 5)
draw(tf_cen_cor_hm)
dev.off()


####


save(
  broadcols,
  nw_names,
  file = paste0("data/rda/ananse_graph_utils.rda")
)

save(
  lg,
  tf_cen,
  tf_cen_rel,
  file = paste0("data/rda/ananse_graph_analysis.rda")
)

########

## ANANSE INFLUENCE STUFF
## tables of influence
infl_dir <- paste0("./scANANSE/analysis_11_2025/influence/")

infl_tables <- 
  paste0(
    infl_dir,
    list.files(path = infl_dir, pattern = ".tsv")
  )

infl_tables= infl_tables[grepl("15.1|45|35",infl_tables)]
infl_tables=infl_tables[c(grepl("_average",infl_tables))]
infl_tables=infl_tables[!c(grepl("_diffnetwork",infl_tables))]


names(infl_tables) <- gsub(infl_dir,"",gsub("anansesnake_", "",gsub(".tsv","",infl_tables)))


###

ref_colors <- c("#368236","#0e6655", "#229954", "#138d75")

lighten_2steps <- function(color, n = 4) {
  lighten_seq <- seq(0, 0.6, length.out = n)  
  lighten(color, amount = lighten_seq)
}

color_list2 <- lapply(ref_colors, lighten_2steps)

all_colors2 <- unlist(color_list2)

broadcols_inf <- all_colors2

list_contrasts=names(infl_tables)
names(broadcols_inf)=list_contrasts
names_order_inf <- names(broadcols_inf)

infls <- setNames(names_order_inf , names_order_inf)


# networks of influence
infl_nws <-
  paste0(
    infl_dir,
    list.files(path = infl_dir, pattern = "_diffnetwork.tsv")
  )

infl_nws= infl_nws[grepl("average",infl_nws)]

names(infl_nws) <- gsub(infl_dir,"",gsub("anansesnake_", "",gsub("_diffnetwork.tsv","",infl_nws)))

infl_nws <- infl_nws[match(names_order_inf, names(infl_nws))]


### load DEGs, : overlay differential expression onto the network nodes
deseq_dir <- paste0("./scANANSE/analysis_11_2025/deseq2/")

DEG_nws <-
  paste0(
    deseq_dir,
    list.files(path = deseq_dir, pattern = ".diffexp.tsv")
  )

DEG_nws= DEG_nws[grepl("15.1|35.|45",DEG_nws)]
DEG_nws= DEG_nws[grepl("average",DEG_nws)]

names(DEG_nws) <- gsub(deseq_dir,"",gsub("danRer11-anansesnake_", "",gsub(".diffexp.tsv","",DEG_nws)))

DEG_nws <- DEG_nws[match(names_order_inf, names(DEG_nws))]

#For the sake of the visualisation, we will be retrieving the top 3 outgoing interactions per TF:

lg_inf <- lapply(1:length(infls), function(x) list())
names(lg_inf) <- infls

for (i in 1:length(lg_inf)){
  message(i)
  inf <- read.table(infl_tables[i], header = TRUE)
  nw <- read.table(infl_nws[i], header = TRUE)
  g <- graph.data.frame(d = nw, directed = TRUE)
  deg=  read.table(DEG_nws[i], header = TRUE)
  
  V(g)$centr_score <- nan_to_zero(closeness(g, mode = "out"))
  V(g)$outdegree <- nan_to_zero(igraph::degree(g, mode = "out"))
  if(length(V(g)$outdegree)>0){
    g_orig <- g
    
    g <- induced_subgraph(g, vids = V(g)$name[V(g)$name %in% inf$factor[1:20]])
    g <- igraph::delete.vertices(g, which(igraph::degree(g) == 0))
    
    comp <-names(lg_inf)[i] #str_split(names(lg_inf)[i],"_")[[1]][1]
    pal <- colorRampPalette(c("#ececec", broadcols_inf[comp]))(5)
    give_col <- circlize::colorRamp2(seq(0,1, len = length(pal)),pal)
    
    V(g)$centr_score <- relativise(V(g)$centr_score)
    V(g)$centr_colour <- give_col(V(g)$centr_score)
    
    V(g)$outdegree <- relativise(V(g)$outdegree)
    V(g)$outdegree_colour <- give_col(V(g)$outdegree)
    
    V(g)$infl_score <- inf$influence_score[match(V(g)$name, inf$factor)]
    V(g)$infl_colour <- give_col(V(g)$infl_score)
    
    # If we want to highlight the D.E. TFs, we can do this:
    diffgenes= deg %>% dplyr::filter(padj<0.05) %>% dplyr::filter(abs(log2FoldChange) > log2(1.5))
    diffgenes.n= rownames(diffgenes)
    if(length(diffgenes.n)>0){
      V(g)$size <- ifelse(V(g)$name %in% diffgenes.n, 2, .6)
    }
    
    
    message(i)
    top_edges <- integer()
    for(tf in V(g)$name[igraph::degree(g,V(g), mode = "out")>0]){
      message(tf)
      e <- top_edges_per_tf(g=g, tf = tf, top = 5)
      top_edges <- c(top_edges,e)
    }
    
    g_ <- subgraph.edges(graph = g, eids = E(g)[top_edges], delete.vertices = TRUE)
    
    lg_inf[[i]]$influence <- inf
    lg_inf[[i]]$diff_g <- g_
    lg_inf[[i]]$diff_g_orig <- g_orig
    
    rm(inf,nw,g,g_,e)
  }
}




## Visualising the influence results

l_inf = lapply(infl_tables,read.delim2)
#l_inf = lapply(l_inf, function(x){})
l_inf = lapply(l_inf, function(x){x$name = translate_ids(x$factor, tfs_all); return(x)})

for(i in 1:16){
  l_inf[[i]]$col = broadcols_inf[[i]]
  l_inf[[i]]$influence_score=as.numeric(l_inf[[i]]$influence_score)
  l_inf[[i]]$factor_fc=as.numeric(l_inf[[i]]$factor_fc)
}

l_inf_plot=l_inf

for(i in names(l_inf_plot)){
  if (nrow(l_inf_plot[[i]])>0) {
    if(grepl("pre",i)== TRUE){
      l_inf_plot[[i]]$factor_fc=as.numeric(as.numeric(l_inf_plot[[i]]$factor_fc)*-1)
      
    }
  }
}

l_inf_plot_comparison=list()
comps=c("15.1_bPAC","15.1_cont","35.0_bPAC","35.0_cont",
        "35.1_bPAC","35.1_cont","45_bPAC","45_cont")


lighten_2steps <- function(color, n = 2) {
  lighten_seq <- seq(0, 0.6, length.out = n)  
  lighten(color, amount = lighten_seq)
}

color_list2 <- lapply(ref_colors, lighten_2steps)

all_colors2 <- unlist(color_list2)

for (i in 0:7) {
  l_inf_plot_comparison[[comps[i+1]]]=rbind(l_inf_plot[[(2*i)+1]],l_inf_plot[[(2*i)+2]])
}

#visualization
pdf(paste0("figures/ananse_influence_scatters_average.pdf"),wi = 9, he = 5)
for (i in 1:8) {
  plot_inf=ggplot(l_inf_plot_comparison[[i]], aes(factor_fc,as.numeric(influence_score))) +
    geom_point(aes(size = direct_targets, colour = as.numeric(influence_score))) +
    xlim(c(-5,5.5))+
    geom_text(
      aes(
        label=ifelse(factor_fc > 0.25|factor_fc < -0.25,as.character(factor),""),
        hjust = 0.5,
        vjust = 2
      ))+
    geom_hline(yintercept = c(0.75), colour = "gray",linetype='dotted')+
    geom_vline(xintercept=c(-log2(1.5),log2(1.5)),linetype='dotted',colour = "gray")+
    theme_classic()+
    labs(title= paste0(comps[i],": postLD vs. preLD"))+
    scale_colour_gradient(low = "#4c4c4c", high = all_colors2[i])
  
  print(plot_inf)
  
}
dev.off()




## Visualising the influence networks

#{r}


pdf(paste0("figures/ananse_diff_networks_symbol_smallsize_average.pdf"), height = 6, width = 6)
for(i in names(lg_inf)){
  set.seed(5678)
  if (nrow(l_inf[[i]])>0) {
    plot(
      main = i,
      lg_inf[[i]]$diff_g,
      layout = layout.graphopt(lg_inf[[i]]$diff_g),
      vertex.color = V(lg_inf[[i]]$diff_g)$outdegree_colour,
      vertex.frame.color = darken(V(lg_inf[[i]]$diff_g)$outdegree_colour,0.5),
      edge.color = rgb(0.1,0.1,0.1,0.25),
      edge.arrow.size = .5,
      vertex.label.family = "Helvetica",
      vertex.label.color = "black",
      vertex.label =translate_ids(V(lg_inf[[i]]$diff_g)$name, tfs_all[,c("gene","Symbol")]),
      vertex.label.cex = 1,
      vertex.size = 8 * V(lg_inf[[i]]$diff_g)$size
    )
  }
}
dev.off()



#We can focus in one example to see that there are TFs that have similarly high scores of influence in different clusters:

#{r}

pdf("./figures/ananse_influence_bPAC_cont_average.pdf", width = 6, height = 6)
p_inf=names(lg_inf)
for (i in 0:7) {
  inf_bPAC_cont <-
    merge(
      lg_inf[[(2*i)+1]]$influence[,c(1,2)],
      lg_inf[[(2*i)+2]]$influence[,c(1,2)],
      by = 1,
      all.x = TRUE,
      all.y = TRUE
    )
  colnames(inf_bPAC_cont) <- c("id","bPAC(dGC-OE)","control")
  rownames(inf_bPAC_cont) <- inf_bPAC_cont$id; inf_bPAC_cont$id <- NULL
  inf_bPAC_cont[is.na(inf_bPAC_cont)] <- 0
  inf_bPAC_cont <- inf_bPAC_cont[rev(order(inf_bPAC_cont$bPAC+inf_bPAC_cont$cont)),]
  inf_bPAC_cont <- inf_bPAC_cont[inf_bPAC_cont$bPAC + inf_bPAC_cont$cont > 0,]
  
  colvec <- setNames(rep("#BDC3C7", nrow(inf_bPAC_cont)),rownames(inf_bPAC_cont))
  
  bPAC_prim=inf_bPAC_cont %>% filter(`bPAC(dGC-OE)` >0.8)%>%filter(control <0.5)
  if(nrow(bPAC_prim)>0){
  for (k in rownames(bPAC_prim)) {
    colvec[k]=all_colors2[i+1]
  }
    dotsizes = setNames(rep(1,nrow(inf_bPAC_cont)), rownames(inf_bPAC_cont))
    dotsizes[names(dotsizes) %in% rownames(bPAC_prim)] <- 2
    
    plot(inf_bPAC_cont, pch = 21, bg = alpha(colvec,.8), col = alpha(darken(colvec,.4),.8), 
         cex = dotsizes, main = paste0(comps[i+1],":influential TFs. postLD and preLD"),ylim=c(0,1))
    text(x = bPAC_prim$`bPAC(dGC-OE)`,
         y = bPAC_prim$control,
         labels = rownames(bPAC_prim),
         pos = 3,
         cex = 0.6
    )
  }else{
    plot(inf_bPAC_cont, pch = 21, bg = alpha(colvec,.8), col = alpha(darken(colvec,.4),.8), 
         cex = dotsizes, main = paste0(comps[i+1],":influence of bPAC and cont."),ylim=c(0,1))
    
  }
  
}
dev.off()


## Co-influential TFs

#{r}
tf_inf_l <- 
  lapply(
    lg_inf,
    function(x){
      x$influence=filter(x$influence, factor_fc > 0) #log2(1.25))
      setNames(x$influence$influence_score,x$influence$factor)
    }
  )
tf_inf <- t(ldply(tf_inf_l,function(s){t(data.frame(unlist(s)))}))
colnames(tf_inf) <- tf_inf[1,]
tf_inf <- as.data.frame(tf_inf[-1,])
tf_inf = data.frame(lapply(tf_inf, function(x) as.numeric(x)),
                    check.names=F, row.names = rownames(tf_inf))
tf_inf[is.na(tf_inf)] = 0

#colnames(tf_inf) <- sapply(str_split(colnames(tf_inf), "_"), function(x) x[1])

tf_inf <- tf_inf[,match(names_order_inf, colnames(tf_inf))]

tf_inf <-
  tf_inf[
    match(
      names(sort(sapply(as.data.frame(t(tf_inf)),function(x){which(x == max(x))[1]}))), #staircase sorting
      rownames(tf_inf)
    ),
  ]

tf_inf <- tf_inf[apply(tf_inf,1,function(x){any(x>.3)}),]
colnames(tf_inf)=sapply(str_split(colnames(tf_inf), "_"), function(x) x[1])

## Defining clusters of co-influential TFs

#{r}
tf_inf_cor <- cor(t(tf_inf))

set.seed(1234)
tf_inf_clu <- hclust(as.dist(1-tf_inf_cor),method = "complete")
plot(tf_inf_clu)
abline(h=c(.5,.7,.8,.9,1,1.5,2), col = divergingx_hcl(7, "Spectral"), lwd = 1.2)


#We tidy and arrange the clustering similar to what we have done in WGCNA previously.

#{r}
cut_thresh <-1.1

tf_inf_clu_hc <- cutree(tf_inf_clu,h=cut_thresh)

broadcols_inf2=broadcols_inf
names(broadcols_inf2)=sapply(str_split(names(broadcols_inf), "_"), function(x) x[1])


coinf_clusters_table <- 
  reorder_modules(
    data.frame(
      tf_inf,
      module = tf_inf_clu_hc
    ),
    order_criterion = names(broadcols_inf2),
    ordering_function = "median",
    thresh_sd = 1 # change for 1.5 if needed
  )

coinf_clusters_table$celltypes=gsub(pattern = "X","",coinf_clusters_table$celltypes)
coinf_clusters_table$col <- 
  broadcols_inf2[match(coinf_clusters_table$celltypes, names(broadcols_inf2))]

names(sort(sapply(as.data.frame(t(tf_inf)),function(x){which(x == max(x))[1]})))

coinf_clusters_table$col[grep(",",coinf_clusters_table$celltypes)] <- 
  sapply(
    coinf_clusters_table$celltypes[grep(",",coinf_clusters_table$celltypes)],
    function(x){
      y = average_cols(broadcols_inf2[unlist(strsplit(x,split=","))])
      return(y)
    }
  )

if(length(coinf_clusters_table$col[is.na(coinf_clusters_table$col)]) !=0){
  coinf_clusters_table$col[is.na(coinf_clusters_table$col)] <- 
    sapply(
      seq(.1,.6, length = length(coinf_clusters_table$col[is.na(coinf_clusters_table$col)])),
      function(x){y=rgb(x,x,x);return(y)}
    )
}

clu_cols <- 
  setNames(
    coinf_clusters_table$col,
    factor(coinf_clusters_table$module_wgcna, levels = unique(coinf_clusters_table$module_wgcna))
  )

tf_inf_cluID <- setNames(
  translate_ids(tf_inf_clu_hc, coinf_clusters_table[,c(2,4)]),
  names(tf_inf_clu_hc)
)

head(coinf_clusters_table)


#Here the overall behaviour of the TFs in these clusters of co-influence

##{r, fig.width = 8, fig.height = 10}
pdf("figures/ananse_influence_tf_clusters_average.pdf", height = 14, width = 8)
par(mfrow=c(5,4))
for(i in coinf_clusters_table$module_wgcna){
  boxplot(
    tf_inf[tf_inf_clu_hc == i,],
    ylim = c(0,1),
    las = 2, xaxt = "n", frame.plot = FALSE,
    main = paste0("TF cluster ", i),
    xlab = "comparison",
    ylab = "influence",
    col = broadcols_inf2, border = darken(broadcols_inf2, .5)
  )
}
par(mfrow = c(1,1))
dev.off()


#We will extract the top coinfluential TFs using correlation to their own module as a proxy:


#{r}
tf_inf_avg <- aggregate(tf_inf, by = list(clu = tf_inf_cluID), FUN = mean)
rownames(tf_inf_avg) <- tf_inf_avg$clu
tf_inf_avg$clu <- NULL

coinf_cor <- 
  sapply(
    rownames(tf_inf),
    function(x){
      y = tf_inf_cluID[x]
      a = as.numeric(tf_inf[x,])
      b = as.numeric(tf_inf_avg[y,])
      z = cor(a,b)
      return(z)
    }
  )

coinf_cor_df <-
  data.frame(
    id = names(coinf_cor),
    cor = coinf_cor,
    clu = tf_inf_cluID
  )

xlsx::write.xlsx(
  coinf_cor_df,
  file = "./outputs/danio_ANANSE_TF_coinfluence_membership_average_average.xlsx",
  sheetName = "tf co-influence cluster membership",
  col.names = TRUE, row.names = FALSE, showNA = TRUE
)

set.seed(1234)
coinf_cor_df %>% 
  group_by(clu) %>% 
  slice_max(order_by = cor, n = 5) %>%
  slice_sample(n=5) -> coinf_cor_df_top

# coinf_kme <- cor(t(tf_inf),t(tf_inf_avg)) # I think we do not use this for anything

tf_inf_top <- tf_inf[coinf_cor_df_top$id,]

head(coinf_cor_df_top)


#Here a heatmap for visualisation:

#  #{r, fig.width = 2.8, fig.height = 8}
# canijo

colnames(tf_inf_top)=sapply(str_split(colnames(tf_inf_top), "_"), function(x) x[1])
tf_inf_hm <-
  Heatmap(
    name = "TF\ninfluence\nscore",
    tf_inf_top,
    col = colorRampPalette(c("#f1f5ff","#b1b8c4","#2c2b46"))(10), # sequential_hcl(10, "Sunset"),
    # show_row_names = FALSE,
    show_column_names = T, column_names_side = "bottom",
    cluster_rows = T,
    clustering_method_rows = "complete",
    clustering_method_columns = "average", #"ward.D2",
    row_split = factor(coinf_cor_df_top$clu, levels = coinf_clusters_table$newname),
    row_labels = translate_ids(rownames(tf_inf_top),tfs_all[,c(1,2)]),
    row_names_gp = gpar(cex = .7),
    row_title_rot = 0,
    show_row_dend = T,
    row_title_gp = gpar(cex=.7),
    top_annotation = 
      HeatmapAnnotation(
        comp = factor(colnames(tf_inf_top), levels = colnames(tf_inf_top)),
        `no. TFs` = anno_barplot(
          sapply(
            lg_inf[match(sapply(str_split(names_order_inf, "_"), function(x) x[1]),
                         sapply(str_split(colnames(tf_inf_top), "_"), function(x) x[1]))],
            function(x){nrow(x$influence)}),
          gp = gpar(fill = broadcols_inf2, border = darken(broadcols_inf2,.5))
        ),
        col = list(comp = setNames(broadcols_inf2,colnames(tf_inf_top))),
        show_legend = FALSE
      ),
    left_annotation = HeatmapAnnotation(
      hc = factor(coinf_cor_df_top$clu, levels = unique(coinf_clusters_table$newname)),
      col = list(
        hc = setNames(clu_cols,coinf_clusters_table$newname[match(names(clu_cols),coinf_clusters_table$module_wgcna)])
      ),
      which = "row",
      show_legend = FALSE
    )
  )



pdf("figures/tf_inf_hm_average.pdf", width =5, height = 20)
draw(tf_inf_hm)
dev.off()


##{r}
colnames(tf_inf)=sapply(str_split(colnames(tf_inf), "_"), function(x) x[1])

tf_inf_hm_supp <-
  Heatmap(
    name = "TF\ninfluence\nscore",
    tf_inf,
    col = colorRampPalette(c("#f1f5ff","#b1b8c4","#2c2b46"))(10), # sequential_hcl(10, "Sunset"),
    show_row_names = FALSE,
    cluster_columns = TRUE,
    show_column_names = T, 
    cluster_rows = tf_inf_clu,
    column_names_side = "bottom",
    
    clustering_method_rows = "complete",#ward.D2",
    top_annotation = 
      HeatmapAnnotation(
        comp = factor(colnames(tf_inf), levels = colnames(tf_inf)),
        `no. TFs` = anno_barplot(
          sapply(
            lg_inf[match(sapply(str_split(names_order_inf, "_"), function(x) x[1]),
                         sapply(str_split(colnames(tf_inf), "_"), function(x) x[1]))],
            function(x){nrow(x$influence)}),
          gp = gpar(fill = broadcols_inf, border = darken(broadcols_inf,.5))
        ),
        col = list(comp = setNames(broadcols_inf,colnames(tf_inf))),
        which = "column"
      ),
    left_annotation =
      HeatmapAnnotation(
        hc = factor(tf_inf_clu_hc, levels = unique(coinf_clusters_table$module_wgcna)),
        col = list(
          hc = clu_cols
        ),
        which = "row",
        show_legend = FALSE
      )
  )
''

pdf("figures/tf_inf_hm_sup_average.pdf",width = 5, height = 20)
draw(tf_inf_hm_supp)
dev.off()

## Visualising co-influential TFs with igraph

##{r, fig.height = 8, fig.width = 8}
g3 <- graph_from_adjacency_matrix(
  cor(t(tf_inf)),mode = "upper",diag = FALSE, weighted = TRUE
)

g3 <- subgraph.edges(g3, eids = E(g3)[E(g3)$weight > .7])
V(g3)$hc_clu <- tf_inf_clu_hc[match(V(g3)$name, names(tf_inf_clu_hc))]
V(g3)$color_hc <- translate_ids(V(g3)$hc_clu, dict = coinf_clusters_table[,c(2,7)])
ran <- .25+relativise(E(g3)$weight)
ran[ran>1] <- 1
E(g3)$color <- rgb(0.1,0.1,0.1,ran)

pdf("figures/g3_average.pdf",width = 8, height = 10)
set.seed(25)
plot(
  g3,
  layout = layout_components(g3),
  vertex.size = 3,
  vertex.label = NA,
  edge.width = E(g3)$weight,
  vertex.color = V(g3)$color_hc,
  vertex.frame.color = darken(V(g3)$color_hc,.4),
  main = "Groups of co-influential TFs"
)
legend(
  'bottomleft',
  legend=coinf_clusters_table$celltypes,
  col= darken(coinf_clusters_table$col, .4),
  pt.bg=coinf_clusters_table$col,pch=21, bty = "n",cex =1
)
dev.off()


## Visualising top target genes of each graph:


#{r}
#load("outputs/rda/ananse_graph_analysis.rda")


list_inf_dfs <- lapply(lg_inf, function(x) x$influence)  # Extract data frames
names(list_inf_dfs)=sapply(str_split(names(list_inf_dfs), "_"), function(x) x[1])
# Create a vector of original sublist names corresponding to each row in data frames
comp <- unlist(lapply(seq_along(list_inf_dfs), function(i) rep(sapply(str_split(names(list_inf_dfs), "_"), function(x) x[1])[i], 
                                                               nrow(list_inf_dfs[[i]]))
))
inf_df_all <- do.call(rbind, list_inf_dfs)
inf_df_all$comp <- comp
inf_df_all$genesymbol <- 
  translate_ids(x = inf_df_all$factor, dict = tfs_all[,c(1,2)])

inf_df_top <- inf_df_all %>% group_by(comp) %>% slice_max(order_by = influence_score, n = 10)


#Here we retrieve the list of top targets per TF of choice, in their respective graphs:

#{r}
list_tgs <- list()

for(i in unique(inf_df_top$comp)){
  
  compa <- unlist(strsplit(i, "_"))[1]
  tfs <- inf_df_top$factor[inf_df_top$comp == i]
  
  list_tgs[[compa]]<- 
    do.call(
      "rbind",
      base::lapply(
        tfs, 
        function(x){
          #if it's in the network
          if(x %in% vertex_attr(lg[[compa]], "name") ){
            y = E(lg[[compa]])[.from(x)]$prob # values
            names(y) = head_of(lg[[compa]], es = E(lg[[compa]])[.from(x)])$name # names
            y = y[y > quantile(y, .95)] # filter
            #if more than 0
            if(length(y)>0){
              z = data.frame( # make a DF out of this
                comp = compa,
                tf = as.character(x),
                tg = names(y),
                value = y,
                value_rel = relativise(y)
              )
              rownames(z) = NULL
              return(z)
            }
          }
        }
      )
    )
}


#Here we retrieve the list of top targets per TF of choice, in their respective graphs:

DF_tgs <- do.call("rbind",list_tgs) # genius

DF_tgs$tf_sym <- translate_ids(x=DF_tgs$tf, dict = tfs_all[,c(1,2)])

DF_tgs$color <-
  broadcols_inf[
    match(
      DF_tgs$comp,
      sapply(str_split(names(broadcols_inf), "_"), function(x) x[1])
    )
  ]

DF_tgs$comp <- factor(DF_tgs$comp, levels = sapply(str_split(names(broadcols_inf), "_"), function(x) x[1]))

DF_tgs <- DF_tgs[order(DF_tgs$comp, DF_tgs$tf),]

DF_tgs$comp_tf <- paste(as.character(DF_tgs$comp), DF_tgs$tf_sym)
DF_tgs$comp_tf <- factor(DF_tgs$comp_tf, levels = unique(DF_tgs$comp_tf))

DF_tgs$name <- DF_tgs$gene
DF_tgs$name[is.na(DF_tgs$name)] <- ""
head(DF_tgs)


library(openxlsx)

write.xlsx(
  DF_tgs,
  file = "./outputs/danio_ANANSE_target_lists_average.xlsx",
  sheetName = "target_lists",
  colNames = T, rowNames = TRUE, showNA = TRUE
)


####TFs for DEG : identify the transcription factors predicted to drive the DEGs
#DE-Tgs
mild_DEGs=list()
for (i in levels(DF_tgs$comp)){
  a = read.table(DEG_nws[paste0(i,"_average")], header = TRUE) 
  a = a %>% dplyr::filter(padj <0.05) %>%dplyr::filter(abs(log2FoldChange)>log2(1.5))
  if(nrow(a)>0){
  mild_DEGs[[i]]=a
}else{ mild_DEGs[[i]]=""}
}


list_DE_tgs <- list()

for(i in unique(inf_df_top$comp)){
  
  compa <- unlist(strsplit(i, "_"))[1]
  DEGs <- rownames(mild_DEGs[[i]])
  
  list_DE_tgs[[compa]]<- 
    do.call(
      "rbind",
      base::lapply(
        DEGs, 
        function(x){
          #if it's in the network
          if(x %in% vertex_attr(lg[[compa]], "name") ){
            y = E(lg[[compa]])[.to(x)]$prob # values
            names(y) = tail_of(lg[[compa]], es = E(lg[[compa]])[.to(x)])$name # names
            y = y[y > quantile(y, .95)] # filter
            #if more than 0
            if(length(y)>0){
              z = data.frame( # make a DF out of this
                comp = compa,
                tg = as.character(x),
                tf = names(y),
                value = y,
                value_rel = relativise(y),
                level_direct=1
              )
              rownames(z) = NULL
              return(z)
            }
          }
        }
      )
    )
}

###add secondary TFs

# target_tf_find(): walk the network to find TFs regulating a target-gene list,
# up to `direct_level` steps upstream (direct + secondary regulators).
target_tf_find=function(detg_list,list_name,direct_level){
  list_name <- list()
  for(i in unique(inf_df_top$comp)){
    
    compa <- unlist(strsplit(i, "_"))[1]
    DEGs <- unique(c(detg_list[[i]]$tg[which(detg_list[[i]]$tg %in% tfs_all$gene)]))
    
    list_name[[compa]]<- 
      do.call(
        "rbind",
        base::lapply(
          DEGs, 
          function(x){
            #if it's in the network
            if(x %in% vertex_attr(lg[[compa]], "name") ){
              y = E(lg[[compa]])[.to(x)]$prob # values
              names(y) = tail_of(lg[[compa]], es = E(lg[[compa]])[.to(x)])$name # names
              y = y[y > quantile(y, .8)] # filter
              #if more than 0
              if(length(y)>0){
                z = data.frame( # make a DF out of this
                  comp = compa,
                  tg = as.character(x),
                  tf = names(y),
                  value = y,
                  value_rel = relativise(y),
                  level_direct=direct_level
                )
                rownames(z) = NULL
                return(z)
              }else{
                print(paste0("no score for", x))
              }
              
            }else{
              print(paste0("no connection for", x))
            }
            
          }
        )
      )
  }
  return(list_name)
}

#2nd tg
list_DE_2tgs=target_tf_find(list_DE_tgs,list_DE_2tgs,2)
#3rd tg
list_DE_3tgs=target_tf_find(list_DE_2tgs,list_DE_3tgs,3)
#4th tg
list_DE_4tgs=target_tf_find(list_DE_3tgs,list_DE_4tgs,4)
#5th tg
list_DE_5tgs=target_tf_find(list_DE_4tgs,list_DE_5tgs,5)

DE_DF_tgs <- do.call("rbind",c(list_DE_tgs,list_DE_2tgs,list_DE_3tgs,list_DE_4tgs,list_DE_5tgs)) # genius

DE_DF_tgs$tf_sym <- translate_ids(x=DE_DF_tgs$tf, dict = tfs_all[,c(1,2)])

DE_DF_tgs$color <-
  broadcols_inf[
    match(
      DE_DF_tgs$comp,
      sapply(str_split(names(broadcols_inf), "_"), function(x) x[1])
    )
  ]

head(DE_DF_tgs)


library(openxlsx)

write.xlsx(
  DE_DF_tgs,
  file = "./outputs/danio_ANANSE_DEG_TF_lists_average.xlsx",
  sheetName = "lists",
  colNames = T, rowNames = TRUE, showNA = TRUE
)

### heatmap TG TF
DE_DF_tgs$value_rel[is.nan(DE_DF_tgs$value_rel)]=0

DE_DF_tgs_mt=list()

for(i in unique(inf_df_top$comp)){
  dt=DE_DF_tgs[which(DE_DF_tgs$comp == i),]
  if(nrow(dt)>0){
  tf_d=unique(dt$tf)
  tg_d=unique(dt$tg)
  
  mt=matrix(0, nrow = length(tf_d),ncol = length(tg_d))
  rownames(mt)=tf_d
  colnames(mt)=tg_d
  
  for (k in 1:nrow(dt)) {
    mt[dt$tf[k],dt$tg[k]]=dt$value_rel[k]
  }
  
  DE_DF_tgs_mt[[i]] <-
    Heatmap(
      name = "TF-TG\nr.score",
      mt,
      col = colorRampPalette(c("#f1f5ff","#b1b8c4","#2c2b46"))(10), # sequential_hcl(10, "Sunset"),
      show_row_names = T,
      show_column_names = T, column_names_side = "bottom",
      cluster_rows = T,
      #clustering_method_rows = "complete",
      #clustering_method_columns = "average", #"ward.D2",
      row_names_gp = gpar(cex = .7),
      #row_title_rot = 0,
      show_row_dend = T,
      row_title_gp = gpar(cex=.7),
      km=3
    )
  
  pdf(paste0("./figures/DE_DF_tgs_mt_hm_",i,"_average.pdf"), width =15, height = 15)
  draw(DE_DF_tgs_mt[[i]])
  dev.off()
  }
}



#This DF can be parsed and used to create lists of plots for each graph separately.

#{r}
ppp <- list()
for(i in levels(DF_tgs$comp)){
  message(i)
  
  d <- DF_tgs[DF_tgs$comp == i, ]
  if(nrow(d)>0){
    p_l <- list()
    for(j in unique(d$tf_sym)){
      d_ = d[d$tf_sym == j,]
      
      l_ <- length(unique(d_$tg[!(d_$tg %in% c(" ","-",""))]))
      
      if(l_ > 30 ){
        set.seed(1234)
        margin_text <- rev(sort(
          sample(unique(d_$tg[!(d_$tg %in% c(" ","-",""))]),30)
        ))  
      } else{
        margin_text <- rev(sort(unique(c(d_$tg[!(d_$tg %in% c(" ","-",""))]))))
      }
      
      margin_text <- margin_text[match(d_$tg[order(d_$value,decreasing = TRUE)], margin_text)]
      margin_text <- margin_text[complete.cases(margin_text)]
      
      d_$in_text <- ifelse(d_$tg %in% margin_text, TRUE, FALSE)
      
      d_ratio_ <- paste0(round(l_ / length(unique(DF_tgs$tg)),2)*100,"%")
      
      set.seed(1234)
      p_l[[j]] <- ggplot(d_ %>% arrange(in_text), aes(x = tf_sym, y = value, fill = color, color = in_text, size = in_text)) +
        geom_jitter(width = 0.3, alpha = 0.7, pch = 21) +
        scale_size_manual(values = c(1.5,2.5))+
        scale_fill_identity() +
        scale_color_manual(values = c("white","black"))+
        labs(#x = "TF Symbol", 
          y = paste0("Score (top 5% targets) (",d_ratio_," in total-TGs.)")) +
        theme_minimal() +
        theme(plot.margin = unit(c(0, 0, 0, 0), "in"))+
        guides(color ="none", size = "none")+
        grid.text(
          label = margin_text, x = .4, y = seq(0.20, 0.95, length.out = length(margin_text)), 
          gp = gpar(fontsize = 6,fontface = "italic", hjust = 1), draw = FALSE
        )
    }
    
    ppp[[i]] <- plot_grid(plotlist = p_l, ncol = 10, nrow = 1)
  }
}

pdf("figures/ananse_top_stripchart_average.pdf",height = 3.5,width = 20)
for(i in names(ppp)){
  print(ppp[[i]])
}
dev.off()


## Saving everything:

#{r}
save(
  lg_inf,
  tf_inf,
  tf_inf_top,
  coinf_clusters_table,
  coinf_cor_df,
  coinf_cor_df_top,
  DF_tgs,
  DE_DF_tgs,
  file = "./data/rda/ananse_coinfluence_average.rda"
)

###### specific influence : TF influence on the DEG programme in the averaged network

## Visualising the influence networks of specific TFs

#In the previous markdown we visualised a subset of the influential networks. 
In this markdown we will visualise the networks of several transcription factors of interest 
that have been previously studied in the literature. 
We have functional knockdown data (bulk RNA-Seq) from these, 
which means that we can dig the validity of our networks by checking the position of the DEGs in the network 
(as having higher interactions etc.).
#We aim at showing the TF of interest and their top interactions with other TFs. 
Color indicates influence on the cell type. 
If the surrounding network of the TF of interest is very dense (more than 50 interactions),
we subset to keep only interactions that are even more top, followed by some high interactions between the neighbour genes only.

#r}
lg_functional <- list()

crema = "crema"
nr4a1 = "nr4a1"
esr2b = "esr2b"
ahr2 = "ahr2"
nr1d2a = "nr1d2a"
ppargc1a="ppargc1a"

tf_fate <-
  list(
    
    crema = c("crema","45.sst1.1.bPAC.post"),
    nr4a1 = c("nr4a1","45.sst1.1.bPAC.post"),
    esr2b = c("esr2b","35.1.avp.bPAC.post"),
    ahr2 = c("ahr2","35.0.avp.crhb.cont.post"),
    nr1d2a = c("nr1d2a","35.0.avp.crhb.bPAC.post"),
    ppargc1a=c("ppargc1a","35.0.avp.crhb.bPAC.post")
    
  )
broadcols_inf2=broadcols_inf
names(broadcols_inf2)=sapply(str_split(names(broadcols_inf2), "_"), function(x) x[1])
names(lg_inf)=sapply(str_split(nameslg_inf, "_"), function(x) x[1])


for(i in names(tf_fate)){
  tf <- tf_fate[[i]][1]
  fate <- tf_fate[[i]][2]
  inf <- lg_inf[[fate]]$influence
  
  g_f <- lg_inf[[fate]]$diff_g_orig
  g_f <- induced_subgraph(g_f, vids = which(V(g_f)$name %in% V(subgraph_from_edges(g_f, E(g_f)[.from(tf)]))$name ))
  g_f <- induced_subgraph(g_f, which(V(g_f)$name %in% danio_all_tf_ananse))
  g_f <- subgraph_from_edges(g_f, which(E(g_f)$weight > 0.5 ))
  
  if( ecount(g_f) > 50 ){
    
    # top 20 from the one of interest
    g2 <- subgraph_from_edges(g_f, eids = E(g_f)[.from(tf)] )
    nn <- min(ecount(g2),50)
    g2 <- subgraph_from_edges(g2, eids = which(E(g2)$weight >= rev(sort(E(g2)$weight))[nn]) )
    
    # some high ones from others?
    g3 <- induced_subgraph(g_f, vids = which(V(g_f)$name %in% V(g2)$name),impl = "create_from_scratch" )
    g3 <- subgraph_with_top_edges_per_tf(g3,top = 3, delete_isolated = TRUE,mode = "both")
    g3 <- delete.vertices(g3, v = V(g3)$name[V(g3)$name %in% tf])
    
    attrs <- rbind(igraph::as_data_frame(g2, "vertices"), igraph::as_data_frame(g3, "vertices")) %>% unique()
    el <- rbind(igraph::as_data_frame(g2), igraph::as_data_frame(g3)) %>% unique()
    new_g <- graph_from_data_frame(el, directed = TRUE, vertices = attrs)
    
    new_g <- induced_subgraph(new_g, vids = which(igraph::degree(new_g,V(new_g), mode = "all") > 1),impl = "create_from_scratch" )
    
    g_f <- new_g
    
    rm(g2,g3,new_g,nn)
    
  } 
  
  pal <- colorRampPalette(c("#ececec", broadcols_inf2[fate]))(5)
  give_col <- circlize::colorRamp2(seq(0,1, len = length(pal)),pal)
  
  V(g_f)$infl_score <- inf$influence_score[match(V(g_f)$name, inf$factor)]
  V(g_f)$infl_colour <- give_col(V(g_f)$infl_score)
  V(g_f)$infl_colour[is.na(V(g_f)$infl_colour)] <- "#ececec"
  
  lg_functional[[i]] <- g_f
  message(i)
}


for(i in names(lg_functional)){
  set.seed(5678)
  plot(
    main = i,
    lg_functional[[i]],
    layout = layout_nicely(lg_functional[[i]]),
    vertex.color = V(lg_functional[[i]])$infl_colour,
    vertex.frame.color = darken(V(lg_functional[[i]])$infl_colour,0.5),
    edge.color = rgb(0.1,0.1,0.1,0.25),
    edge.arrow.size = .5,
    vertex.label.family = "Helvetica",
    vertex.label.color = "black",
    vertex.label = translate_ids(V(lg_functional[[i]])$name, tfs_all[,c("gene","Symbol")]),
    vertex.label.cex = .5,
    vertex.size = 12 #* V(lg_functional[[i]])$size
  )
}


#r}
pdf("./figures/ananse_networks_tfs_functional_summarised_average.pdf", height = 6, width = 6)
for(i in names(lg_functional)){
  set.seed(5678)
  plot(
    main = i,
    lg_functional[[i]],
    layout = layout_nicely(lg_functional[[i]]),
    vertex.color = V(lg_functional[[i]])$infl_colour,
    vertex.frame.color = darken(V(lg_functional[[i]])$infl_colour,0.5),
    edge.color = rgb(0.1,0.1,0.1,0.25),
    edge.arrow.size = .5,
    vertex.label.family = "Helvetica",
    vertex.label.color = "black",
    vertex.label = translate_ids(V(lg_functional[[i]])$name, tfs_all[,c("gene","Symbol")]),
    vertex.label.cex = .5,
    vertex.size = 12 #* V(lg_functional[[i]])$size
  )
}
dev.off()


#Please refer to the re-analysis markdowns to see the associated boxplots of interactions between these TFs and the down-regulated genes detected in those previous works in the literature.

## Showing networks of co-influential TFs

#We are going to visualise now the network of HNF4 and alx3, two transcription factors that were detected as highly co-influential of two major groups of cell types.

#r}
hnf4 <- "h1SMcG0019688"
alx3 <- "h1SMcG0021040"

lg_coinf <- list()

coinf_fate <-
  list(
    h_pha = c(hnf4,"phagocytes_neoblast"),
    h_par = c(hnf4,"parenchyma_neoblast"),
    h_bg = c(hnf4,"basalgoblet_neoblast"),
    a_neu = c(alx3,"neuron_neoblast"),
    a_mus = c(alx3,"muscle_neoblast"),
    a_sec = c(alx3,"secretory_neoblast")
  )



#Here we create these graphs following exactly the same criterion as before.

#r}
for(i in names(coinf_fate)){
  tf <- coinf_fate[[i]][1]
  fate <- coinf_fate[[i]][2]
  inf <- lg_inf[[fate]]$influence
  
  g_f <- lg_inf[[fate]]$diff_g_orig
  g_f <- induced_subgraph(g_f, vids = which(V(g_f)$name %in% V(subgraph_from_edges(g_f, E(g_f)[.from(tf)]))$name ))
  g_f <- induced_subgraph(g_f, which(V(g_f)$name %in% smed_tfs$id))
  g_f <- subgraph_from_edges(g_f, which(E(g_f)$weight > 0.5 ))
  
  if( ecount(g_f) > 30 ){
    
    # v_add <- E()
    # top 20 from the one of interest
    g2 <- subgraph_from_edges(g_f, eids = E(g_f)[.from(tf)] )
    nn <- min(ecount(g2),30)
    g2 <- subgraph_from_edges(g2, eids = which(E(g2)$weight >= rev(sort(E(g2)$weight))[nn]) )
    
    # some high ones from others?
    g3 <- induced_subgraph(g_f, vids = which(V(g_f)$name %in% V(g2)$name),impl = "create_from_scratch" )
    g3 <- subgraph_with_top_edges_per_tf(g3,top = 3, delete_isolated = TRUE,mode = "both")
    g3 <- delete.vertices(g3, v = V(g3)$name[V(g3)$name %in% tf])
    
    attrs <- rbind(as_data_frame(g2, "vertices"), as_data_frame(g3, "vertices")) %>% unique()
    el <- rbind(as_data_frame(g2), as_data_frame(g3)) %>% unique()
    new_g <- graph_from_data_frame(el, directed = TRUE, vertices = attrs)
    
    g_f <- new_g
    
    rm(g2,g3,new_g,nn)
    
  } 
  
  pal <- colorRampPalette(c("#ececec", broadcols_inf[gsub("_neoblast","",fate)]))(5)
  give_col <- circlize::colorRamp2(seq(0,1, len = length(pal)),pal)
  
  V(g_f)$infl_score <- inf$influence_score[match(V(g_f)$name, inf$factor)]
  V(g_f)$infl_colour <- give_col(V(g_f)$infl_score)
  V(g_f)$infl_colour[is.na(V(g_f)$infl_colour)] <- "#ececec"
  
  lg_coinf[[i]] <- g_f
  message(i)
}


Here's the plot:
  
  #r, fig.width = 6, fig.height = 6}
for(i in names(lg_coinf)){
  set.seed(5678)
  plot(
    main = i,
    lg_coinf[[i]],
    layout = layout.graphopt(lg_coinf[[i]]),
    vertex.color = V(lg_coinf[[i]])$infl_colour,
    vertex.frame.color = darken(V(lg_coinf[[i]])$infl_colour,0.5),
    edge.color = rgb(0.1,0.1,0.1,0.25),
    edge.arrow.size = .5,
    vertex.label.family = "Helvetica",
    vertex.label.color = "black",
    vertex.label = translate_ids(V(lg_coinf[[i]])$name, tfs_all[,c("gene","Symbol")]),
    vertex.label.cex = .5,
    vertex.size = 8 #* V(lg_coinf[[i]])$size
  )
}


#r}
pdf("graphics/ananse_networks_hnf_alx_average.pdf", wi = 5, he = 5)
for(i in names(lg_coinf)){
  set.seed(5678)
  plot(
    main = i,
    lg_coinf[[i]],
    layout = layout.graphopt(lg_coinf[[i]]),
    vertex.color = V(lg_coinf[[i]])$infl_colour,
    vertex.frame.color = darken(V(lg_coinf[[i]])$infl_colour,0.5),
    edge.color = rgb(0.1,0.1,0.1,0.25),
    edge.arrow.size = .5,
    vertex.label.family = "Helvetica",
    vertex.label.color = "black",
    vertex.label = translate_ids(V(lg_coinf[[i]])$name, tfs_all[,c("gene","Symbol")]),
    vertex.label.cex = .5,
    vertex.size = 8 #* V(lg_coinf[[i]])$size
  )
}
dev.off()

