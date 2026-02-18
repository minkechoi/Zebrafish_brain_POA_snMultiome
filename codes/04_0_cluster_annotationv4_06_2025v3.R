
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc()
# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(SCP)
library(plyr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(harmony)
library(colorspace)
library(viridis)
library(igraph)
library(Signac)
library(GenomicRanges)
library(rtracklayer)
library(ensembldb)
library(EnrichedHeatmap)
library(ComplexHeatmap)


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



# rds file loading from previous scripts
ss= readRDS(file=paste0("./data/rds/step3_norm_harmony_wnn_RNA_ft_macs",r.variable,".rds"))
#ss=ft_comb_seurat


####selected dim and clustering = res 
#set cell identity, based on the clustering results, res
DefaultAssay(ss)="SCT"
ss$seurat_clusters=factor(ss$seurat_clusters, str_sort(unique(ss$seurat_clusters,numeric = T)))
Idents(ss)=ss$seurat_clusters

# cluster color setting ---------------------------------------------------
library(magrittr)
library(RColorBrewer)

cl_colors <- 
  c(divergingx_hcl(8,"ArmyRose"),
    divergingx_hcl(11,"RdYlBu"),
    divergingx_hcl(7,"Zissou 1"),
    divergingx_hcl(11,"Spectral"),
    divergingx_hcl(8,"Fall"),
    sequential_hcl(8,"Hawaii"),
    divergingx_hcl(8,"Cividis")
  )

num_clusters <- length(unique(ss$seurat_clusters))

set.seed(368)
cols <- sample(unname(cl_colors),num_clusters)

# cell counts for each cluster --------------------------------------------
##cell count, codes from "Alberto Perez-Posada @apposada"
df <- 
  data.frame(
    table(
      ss$orig.ident, ss$seurat_clusters
    )
  )

# Rename the columns of the data frame
colnames(df) <- 
  c("Library", "Cluster", "Count")

# Normalize the counts in each cluster
df$logCount <- 
  log(df$Count,10) / log(tapply(df$Count, df$Cluster, sum)[df$Cluster],10)

df$CountNorm <- df$Count / tapply(df$Count, df$Cluster, sum)[df$Cluster]

##visulization
# Create a stacked barplot raw number
ncells_bar <- ggplot(df, aes(x = Cluster, y = Count, fill = Library)) +
  geom_col(position = "stack") +
  scale_fill_brewer(
    type = "div",
    palette = "Spectral",
    direction = 1,
    aesthetics = "fill"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) +
  ylab("Number of cells") +
  guides(fill = FALSE)

# Create a stacked barplot log number
logncells_bar <- ggplot(df, aes(x = Cluster, y = logCount, fill = Library)) +
  geom_col(position = "stack") +
  scale_fill_brewer(
    type = "div",
    palette = "Spectral",
    direction = 1,
    aesthetics = "fill"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)) +
  ylab("log(Number of cells)") +
  guides(fill = FALSE)

# Create a stacked barplot norm number
normncells_bar <- ggplot(df, aes(x = Cluster, y = CountNorm, fill = Library)) +
  geom_col(position = "stack") +
  scale_fill_brewer(
    type = "div",
    palette = "Spectral",
    direction = 1,
    aesthetics = "fill"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ylab("Proportion of cells")

# Create a grid of the aligned plots
grid <- plot_grid(
  ncells_bar, logncells_bar, normncells_bar,
  nrow = 3,
  align = "v",
  axis = "tb",
  labels = c("A", "B", "C")
)


# Plot ncells per cluster
dir.create(paste0("./figures/cell_type"))

tiff(paste0("./figures/cell_type/var",r.variable,"_cell_numbers_ratio.tiff"),
     width = 24,height = 30,units = "cm", res = 300,compression = "lzw",bg = NA)
print(grid)

dev.off()
# Finding markers ---------------------------------------------------------


ss=PrepSCTFindMarkers(ss)
#identifiy markers
all.markers= FindAllMarkers(
  ss, assay = "SCT",
  group.by = "seurat_clusters",
  logfc.threshold = 0.1,
  test.use = "wilcox",#"DESeq2","MAST","wilcox"
  #slot = "data",
  min.pct = 0.01,
  min.diff.pct = -Inf,
  node = NULL,
  verbose = TRUE,
  only.pos = T,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  mean.fxn = NULL,
  fc.name = NULL,
  base = 2,
  return.thresh = 0.01,
  densify = FALSE
)
#top50s
top100_mk= all.markers %>% 
  dplyr::filter(abs(avg_log2FC)>1) %>%
  dplyr::filter(p_val_adj <0.05) %>%
  group_by(cluster) %>%
  top_n(100, wt=abs(avg_log2FC))

dir.create(paste0("./outputs/cell_type"))

write.csv(all.markers,paste0("./outputs/cell_type/all.markers_var",r.variable,".csv"))


#Identifying target cells of interest ------------------------------------------------------


## neuropeptidergic cells
#clusters for np+ cells over 100

neuropep=c("oxt","avp","th","th2","sst1.1","galn","npvf","fshb","agrp","pmch","pomc",
           "hcrt","npffl","nmbb","nmba","edn1","edn2","edn3b","calca","prlh2","kiss1",
           "npy","pyya","vip","ccka","penka","penkb","nts","crhb","trh","tshba",
           "tshbb","gnrh3","gnrh2","lhb","cga")

top_nps= intersect(neuropep,top100_mk$gene)
top_nps_clusters=unique(top100_mk$cluster[which(top100_mk$gene %in% top_nps)])

#neuropeptide cell clusters
np_clusters=c()

for (i in top_nps_clusters) {
  n_cell=length(ss$seurat_clusters[which(ss$seurat_clusters ==i)])
  if (n_cell > 50) {
    top_nps_over50 = intersect(top_nps,top100_mk$gene[which(top100_mk$cluster ==i)])
    mkt=top100_mk %>% 
      dplyr::filter(cluster == i) %>% 
      dplyr::filter(gene %in% top_nps_over50)%>% 
      dplyr::filter(pct.1 > 0.1)
    
    if (nrow(mkt)>0) {
      print(paste0("seurat_clusters:cluster_",i," = possible cluster for ", mkt$gene ))
      mkt=mkt[which((mkt$pct.1 *n_cell) >50),]
      np_clusters=c(np_clusters,paste(mkt$gene,mkt$pct.1,mkt$cluster,(n_cell*mkt$pct.1),sep = "_"))
      #np_clusters=rbind(np_clusters,mkt)
    }
  }
}

np50_clusters_sp=str_split(np_clusters,pattern = "_")
np50_clusters_tb = data.frame(matrix(unlist(np50_clusters_sp), nrow=length(np50_clusters_sp), byrow=TRUE))
colnames(np50_clusters_tb)=c("gene","pct.1","cluster","ncell")
np50_clusters_tb=dplyr::arrange(np50_clusters_tb,desc(ncell))
npover50_in_cl=unique(np50_clusters_tb$gene)


# Create vector with required format
result <- sapply(split(np50_clusters_tb, np50_clusters_tb$cluster), function(cluster) {
  genes <- unique(cluster$gene)
  if (length(genes) == 1) {
    paste0(cluster$cluster[1], "_", genes)
  } else {
    paste0(cluster$cluster[1], "_", paste(genes, collapse = "/"))
  }
})
np50_clusters <- as.vector(result)

np50_clusters_mask=as.character(ss[["seurat_clusters"]][,1])
spc_np50_np_clusters=str_split(np50_clusters,"_")
for (i in 1:length(np50_clusters)) {
  
  np50_clusters_mask[np50_clusters_mask == spc_np50_np_clusters[[i]][1] ] = paste0(np50_clusters[i])
  
  print(i)
}
ss[["np50_clusters_mask"]]=np50_clusters_mask

###visualization
np50_clusters_mask_simple=np50_clusters_mask
np50_clusters_mask_simple[!grepl("_",np50_clusters_mask_simple)]="."
ss[["np50_clusters_mask_simple"]]=np50_clusters_mask_simple
ss$np50_clusters_mask_simple=factor(np50_clusters_mask_simple,levels = sort(unique(np50_clusters_mask_simple)))

p6=CellDimPlot(
  srt = ss, group.by = "np50_clusters_mask_simple",
  seed = 0,pt.size = 0.7,label.size = 8,
  #palcolor = np_col,
  reduction = umap, theme_use = "theme_blank",#legend.position = "none",
  label = T,label_insitu = T,label_repel = T
)

np_col=p6[["plot_env"]][["colors"]]
np_col["."]="lightgray"
p6=DimPlot(ss,group.by = "np50_clusters_mask_simple", seed = 0,
           reduction = umap,cols = np_col,alpha = 0.8,pt.size = 0.7,
           label = F, label.size =3,stroke.size = 0.5,
           order=T
)+ ggtitle(NULL)& NoLegend()

p6=LabelClusters(p6, id = "np50_clusters_mask_simple", repel = T,  
                 fontface = "bold", color = "black",
                 box = F #,label.padding=0.5,label.size=1,label.r=0
                #,position = "nearest"
)


pdf(paste0("./figures/cell_type/var",r.variable,"_np50_mask_umap.pdf"),
     width = 5,height = 5)
print(p6)

dev.off()


#ncell > 100
np100_clusters=c()

for (i in top_nps_clusters) {
  n_cell=length(ss$seurat_clusters[which(ss$seurat_clusters ==i)])
  if (n_cell > 100) {
    top_nps_over100 = intersect(top_nps,top100_mk$gene[which(top100_mk$cluster ==i)])
    mkt=top100_mk %>% 
      dplyr::filter(cluster == i) %>% 
      dplyr::filter(gene %in% top_nps_over100)%>% 
      dplyr::filter(pct.1 > 0.1)
    
    if (nrow(mkt)>0) {
      print(paste0("seurat_clusters:cluster_",i," = possible cluster for ", mkt$gene ))
      mkt=mkt[which((mkt$pct.1 *n_cell) >100),]
      np100_clusters=c(np100_clusters,paste(mkt$gene,mkt$pct.1,mkt$cluster,(n_cell*mkt$pct.1),sep = "_"))
      #np100_clusters=rbind(np100_clusters,mkt)
    }
  }
}

np100_clusters_sp=str_split(np100_clusters,pattern = "_")
np100_clusters_tb = data.frame(matrix(unlist(np100_clusters_sp), nrow=length(np100_clusters_sp), byrow=TRUE))
colnames(np100_clusters_tb)=c("gene","pct.1","cluster","ncell")
np100_clusters_tb=dplyr::arrange(np100_clusters_tb,desc(ncell))
npover100_in_cl=unique(np100_clusters_tb$gene)

pdf(paste0("./figures/cell_type/var",r.variable,"_np_FeaturePlot.pdf"),
     width = (10*round(length(npover100_in_cl)/2)),
     height = 15*round(length(npover100_in_cl)/(round(length(npover100_in_cl)/2)))
)
p=FeaturePlot(ss,
               features = npover100_in_cl,
               pt.size = 0.7, cols = c("grey80","red"),
               reduction=paste0(umap),
               order = T,min.cutoff = 0.8,
               ncol=round(length(npover100_in_cl)/2)) & NoAxes()& NoLegend()

print(p)
dev.off()

# Create vector with required format
result <- sapply(split(np100_clusters_tb, np100_clusters_tb$cluster), function(cluster) {
  genes <- unique(cluster$gene)
  if (length(genes) == 1) {
    paste0(cluster$cluster[1], "_", genes)
  } else {
    paste0(cluster$cluster[1], "_", paste(genes, collapse = "/"))
  }
})

# Convert to a vector
np100_np_clusters <- as.vector(result)


np100_np_clusters_mask=as.character(ss[["seurat_clusters"]][,1])
spc_np100_np_clusters=str_split(np100_np_clusters,"_")
for (i in 1:length(np100_np_clusters)) {
  
np100_np_clusters_mask[np100_np_clusters_mask == spc_np100_np_clusters[[i]][1] ] = paste0(np100_np_clusters[i])

  print(i)
}

ss[["np100_clusters_mask"]]=np100_np_clusters_mask

###visualization
np100_clusters_mask_simple=np100_np_clusters_mask
np100_clusters_mask_simple[!grepl("_",np100_clusters_mask_simple)]="_"
ss[["np100_clusters_mask_simple"]]=np100_clusters_mask_simple
ss$np100_clusters_mask_simple=factor(np100_clusters_mask_simple,levels = sort(unique(np100_clusters_mask_simple)))

p5=CellDimPlot(
  srt = ss, group.by = "np100_clusters_mask_simple", split.by = "orig.ident", 
  seed = 0,pt.size = 0.7,label.size = 8,
  #palcolor = np_col,
  reduction = umap, theme_use = "theme_blank",legend.position = "none",
  label = F,label_insitu = T,label_repel = T
)

p6=CellDimPlot(
  srt = ss, group.by = "np100_clusters_mask_simple",
  seed = 0,pt.size = 0.7,label.size = 8,
  #palcolor = np_col,
  reduction = umap, theme_use = "theme_blank",legend.position = "none",
  label = T,label_insitu = T,label_repel = T
)

np_col=p6[["plot_env"]][["colors"]]



pdf(paste0("./figures/cell_type/var",r.variable,"_np_mask_umap.pdf"),
     width = 30,height = 30)
 print(p6)

dev.off()

pdf(paste0("./figures/cell_type/var",r.variable,"_np_mask_umap_sp.pdf"),
     width = 30,height = 30)
print(p5)

dev.off()

##### subclustering selection

ft.all.marker=all.markers %>% 
  dplyr::filter(pct.1>0.4)%>% 
  dplyr::filter(avg_log2FC >1)%>% 
  dplyr::filter(p_val_adj <0.01)

low.mk.cl=table(ft.all.marker$cluster)
low.mk.cl=names(low.mk.cl[low.mk.cl<5])

#Find cluster need sub clustering
#selected_cls=subcl_rq(ss)

top_npc=unique(np100_clusters_tb$cluster)

selected_cls=str_sort(unique(c(low.mk.cl,top_npc)),numeric = T)
#selected_cls
#[1] "3"  "6"  "9"  "11" "14" "15" "24" "35" "38" "40" "44" "45" "48" "49"

reso=c(0.02,0.05,0.1,0.2,0.35,0.5)
Idents(ss)=ss$seurat_clusters
dir.create(paste0("./figures/cell_type/sub_clustering"))

for (i in selected_cls) {
  s.cluster=i
  for (j in reso) {
    
    graphname= "wsnn" #"RNA_snn" "SCT_nn"
    subc.name=paste0("sub_",s.cluster)
    
    ss=FindSubCluster(
      ss,
      cluster=s.cluster,
      graph.name=graphname,
      subcluster.name = subc.name,
      resolution = j,
      algorithm =2
    )
    
    sst= subset(ss, seurat_clusters == s.cluster)
    Idents(sst)=sst[[subc.name]][,1]
    a=length(unique(Idents(sst)))
    p=DimPlot(sst, seed = 0,
              reduction=umap, label=T) & NoAxes()
    tiff(paste("./figures/cell_type/sub_clustering/umap",r.variable,"sub_cl",
               i,"res",j,"ncl",a,
               ".tiff",
               sep = "_"),
         width = 13,height = 10,units = "cm", res = 300,compression = "lzw",bg = NA)
    print(p)
    dev.off()
    sst=""
    gc()
  }
}

##### check the resolution that you want before moving the next step
#selected_cls
#[1] "3"  "6"  "9"  "11" "14" "15" "24" "35" "38" "40" "44" "45" "48" "49"

re_selected_cls=c(6,14,15,24,35,40)
re_selected_res=c(0.1,0.1,0.1,0.1,0.1,0.35)

for (i in 1:length(re_selected_cls)) {
  graphname= "wsnn" #"RNA_snn" "SCT_nn"
  s.cluster=re_selected_cls[i]
  subc.name=paste0("sub_",s.cluster)
  
  ss=FindSubCluster(
    ss,
    cluster=s.cluster,
    graph.name=graphname,
    subcluster.name = subc.name,
    resolution = re_selected_res[i],
    algorithm =2
  )
  
}

#manual correction and rename of subcluster and neuroendo cluster
#selected_cls
#[1] "3"  "6"  "9"  "11" "14" "15" "24" "35" "38" "40" "44" "45" "48" "49"

re_selected_cls
#[1]  6 14 15 24 35 40

#6: galn/crhb = c("6_0"), "6_1"
#14: galn =c("14_0"), prlh2 = c("14_1")
#15: trh =c("15_1"), "15_0"
#24: sst1.1/nts = c("24_0"), "24_1","24_2" 
#35: avp.crhb = c("35_0"), avp = c("35_1") 
#40: th= c("40_0","40_1"), "40_2"



#6: galn/crhb = c("6_0"), "6_1"
i=6
re=as.character(ss[[paste0("sub_",i)]][,1])
re[re %in% c("6_0")]="6.0_galn.crhb"
re[re %in% c("6_1")]="6.1_"
ss[[paste0("sub_",i)]][,1]=re


#14: galn =c("14_0"), prlh2 = c("14_1")
i=14
re=as.character(ss[[paste0("sub_",i)]][,1])
re[re %in% c("14_0")]="14.0_galn"
re[re %in% c("14_1")]="14.1_prlh2"
ss[[paste0("sub_",i)]][,1]=re

#15: trh =c("15_1", "15_0")
i=15
re=as.character(ss[[paste0("sub_",i)]][,1])
re[re %in% c("15_1")]="15.1_trh"
re[re %in% c("15_0")]="15.0_trh"
ss[[paste0("sub_",i)]][,1]=re

#24: sst1.1/nts = c("24_0"), "24_1","24_2" 
i=24
re=as.character(ss[[paste0("sub_",i)]][,1])
re[re %in% c("24_0")]="24.0_sst1.1.nts"
re[re %in% c("24_1")]="24.1_"
re[re %in% c("24_2")]="24.2_"
ss[[paste0("sub_",i)]][,1]=re

#35: avp.crhb = c("35_0"), avp = c("35_1") 
i=35
re=as.character(ss[[paste0("sub_",i)]][,1])
re[re %in% c("35_0")]="35.0_avp.crhb"
re[re %in% c("35_1")]="35.1_avp"
ss[[paste0("sub_",i)]][,1]=re

#40: th= c("40_0","40_1"), "40_2"
i=40
re=as.character(ss[[paste0("sub_",i)]][,1])
re[re %in% c("40_0")]="40.0_th"
re[re %in% c("40_1")]="40.1_th"
re[re %in% c("40_2")]="40.2_"
ss[[paste0("sub_",i)]][,1]=re



### reploting subsetts 
for (i in 1:length(re_selected_cls)) {
  sst= subset(ss,seurat_clusters == re_selected_cls[i])
  subc.name=paste0("sub_",re_selected_cls[i])
  a=length(unique(Idents(sst)))
  p=DimPlot(sst, seed = 0,group.by = subc.name,
            reduction=umap, label=T) & NoAxes() & NoLegend()
  tiff(paste0("./figures/cell_type/re_",r.variable,"sub_cl",
              re_selected_cls[i],r.variable,re_selected_res[i],"alg1","ncl",a,
              ".tiff"),
       width = 10,height = 10,units = "cm", res = 300,compression = "lzw",bg = NA)
  print(p)
  sst=NULL
  dev.off()
}


#merging subclusters

metaset=ss[["seurat_clusters"]]

for (i in 1:length(re_selected_cls)) {
  subc.name=paste0("sub_",re_selected_cls[i])
  mset=data.frame(subc.name=ss[[subc.name]])
  metaset=cbind(metaset,mset)
}

# Function to find the unique element(s) for each row
get_unique_or_common <- function(row) {
  unique_vals <- unique(row)
  
  if (length(unique_vals) == 1) {
    return(as.character(unique_vals))  # single value
  } else {
    vals_with_underscore <- unique_vals[grepl("_", unique_vals)]
    return(paste(vals_with_underscore, collapse = ";"))  # collapse to one string
  }
}

# Apply function row-wise
metaset$merged_sub <- apply(metaset[, -1], 1, get_unique_or_common)

ss[["merged_sub"]]=metaset$merged_sub


#kiss1 3,9, 48: th2, 41: nppal
re=as.character(ss[["merged_sub"]][,1])
re[re %in% c("3")]="3_kiss1"
re[re %in% c("9")]="9_kiss1"
re[re %in% c("44")]="44_oxt"
re[re %in% c("45")]="45_sst1.1"
re[re %in% c("48")]="48_th2.th"
re[re %in% c("49")]="49_npy.sst"
re[re %in% c("41")]="41_nppal"

ss[["merged_sub"]]=re


####additional cell number filter (celll min. > 50 per cluster)

lowncells=names(table(ss$merged_sub)[which(table(ss$merged_sub) <100)])
if (length(lowncells)>0) {
  ss=subset(ss, merged_sub %!in% lowncells)
}


#all_marker_merged sub
#identifiy markers
ss=PrepSCTFindMarkers(ss)

all.merged.markers= FindAllMarkers(
  ss, assay = "SCT",group.by = "merged_sub",
  logfc.threshold = 0.1,
  test.use = "wilcox",#"DESeq2","MAST","wilcox"
  #slot = "data",
  min.pct = 0.01,
  min.diff.pct = -Inf,
  node = NULL,
  verbose = TRUE,
  only.pos = T,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  mean.fxn = NULL,
  fc.name = NULL,
  base = 2,
  return.thresh = 0.01,
  densify = FALSE
)
#top50s
top100.merged_mk= all.merged.markers %>% 
  dplyr::filter(abs(avg_log2FC)>1) %>%
  dplyr::filter(p_val_adj <0.05) %>%
  group_by(cluster) %>%
  top_n(100, wt=abs(avg_log2FC))

write.csv(all.merged.markers,paste0("./outputs/cell_type/all.merged.markers_",
                                    "var",r.variable,
                                    ".csv"))



# cluster annotation form refs --------------------------------------------
#ref_cell: markerset1_scheir_lab (Shafer et al., 2022, Nat Ecol Evol.)
po_marker_genes=vroom::vroom("../refs/Supplemental_data/3-marker_gene_lists/Drerio_zebrafish.markers.csv")
po_marker_genes=tibble::column_to_rownames(po_marker_genes, var = "...1")
po_marker_genes[,"gene"]=sapply(strsplit(as.character(row.names(po_marker_genes)), "\\.\\.\\."), `[`, 1)

po_sub_marker_genes=vroom::vroom("../refs/Supplemental_data/3-marker_gene_lists/Drerio_zebrafish.markers.sub.csv")
po_sub_marker_genes=tibble::column_to_rownames(po_sub_marker_genes, var = "...1")
po_sub_marker_genes[,"gene"]=sapply(strsplit(as.character(row.names(po_sub_marker_genes)), "\\.\\.\\."), `[`, 1)

sig.markers_PO= po_marker_genes %>% dplyr::filter(p_val_adj <0.05) %>% dplyr::filter(abs(avg_logFC) > 1)
sig.markers_PO_sub= po_sub_marker_genes %>% dplyr::filter(p_val_adj <0.05) %>% dplyr::filter(abs(avg_logFC) > 1)

#ref_cell2: markerset2_PanglaoDB
PanglaoDB= vroom::vroom("./data/PanglaoDB_markers_27_Mar_2020.tsv")
PanglaoDB_Br=PanglaoDB %>% 
  dplyr::filter(`canonical marker`== 1) %>%
  dplyr::filter(organ == "Brain") %>%
  dplyr::filter(`ubiquitousness index` < 0.05)

db_br_celltype=unique(PanglaoDB_Br$`cell type`)

zfin_human_todanio = vroom::vroom(file = "../refs/annotations/ZFIN/human_orthos_2024.04.24.txt",
                                  delim = "\t", 
                                  col_names = c("ZFIN_ID","ZFIN_Symbol","ZFIN_Name","Human_Symbol",
                                                "Human_Name","OMIM_ID","Gene_ID","HGNC_ID","Evidence","Pub_ID"))

PanglaoDB_zfin_gene = unique(zfin_human_todanio$ZFIN_Symbol[c(which(zfin_human_todanio$Human_Symbol %in% PanglaoDB_Br$`official gene symbol`))])
join_key=join_by("official gene symbol"=="Human_Symbol" )
zfin_PanglaoDB_Br=left_join(PanglaoDB_Br,zfin_human_todanio,by = join_key)
zfin_PanglaoDB_Br_ft=zfin_PanglaoDB_Br[,c(3,16)]
zfin_PanglaoDB_Br_ft=unique(na.omit(zfin_PanglaoDB_Br_ft))
zfin_PanglaoDB_Br_ft$`cell type`[zfin_PanglaoDB_Br_ft$`cell type` =="Immature neurons"]= "IN"
zfin_PanglaoDB_Br_ft$`cell type`[zfin_PanglaoDB_Br_ft$`cell type` =="Neural stem/precursor cells"]= "NS"
zfin_PanglaoDB_Br_ft$`cell type`[zfin_PanglaoDB_Br_ft$`cell type` =="Radial glia cells"]= "RG"
zfin_PanglaoDB_Br_ft$`cell type`[zfin_PanglaoDB_Br_ft$`cell type` =="Oligodendrocyte progenitor cells"]= "OPCs"


# cluster_enrichment_test -------------------------------------------------

### cell enrichment test, Universal enrichment analysis
# https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html

library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(org.Dr.eg.db)
library(msigdbr)

#significant mark genes
sig_markers04=all.merged.markers %>% 
  dplyr::filter(abs(avg_log2FC ) > log2(1.5)) %>% 
  dplyr::arrange(desc(avg_log2FC ))

table(sig_markers04$cluster)

### ref_genelist and cell markers
ref_cell1 = data.frame("cluster"=gsub("\\/","_",sig.markers_PO$cluster),
                       "gene"=sig.markers_PO$gene)
ref_cell2 = zfin_PanglaoDB_Br_ft

rc_list=c("ref_cell1","ref_cell2")

for (rc in rc_list) {
  ref_cell=get(rc)

  upenrich_result=list()
  notlisted=list()
  for (i in unique(sig_markers04$cluster)) {
    print(i)
    mkr = sig_markers04 %>% dplyr::filter(avg_log2FC > log2(1.5)) %>% 
      dplyr::filter(cluster == i)%>% dplyr::arrange(desc(avg_log2FC))
    if(nrow(mkr)>=10){
      upenrich= enricher(mkr$gene, TERM2GENE = ref_cell)
      if (nrow(upenrich@result)>=1) {
        upenrich@result$in.cluster= i
        upenrich_result[[i]]=upenrich@result
      }else{
      print(paste("not enough significant marker for cluster",i))
      notlisted[[i]]=i
    }
    }
  }
  tname=paste0(rc,"_upenrich_result")
  upenrich_result=do.call("rbind", upenrich_result)
  assign(tname,upenrich_result)
  #write.csv(upenrich_result, paste0("./outputs/cell_type/enrichr/up_enrich_cell_clusters",rc,"_.csv" ))
  
}

####annotation voting##
#####
###annotation voting code modifed from "Alberto Perez-Posada @apposada"
#with ref1
#rc="ref_cell1"
#rc="ref_cell2"
for (rc in rc_list) {
  
  tname=paste0(rc,"_upenrich_result")
  df=get(tname)
  df=df%>% dplyr::filter(qvalue <0.05) %>% dplyr::arrange(qvalue)
  
  classifications <- 
    data.frame(
      cluster = df$in.cluster,
      predicted = rownames(df)
    )
  
  classifications_table <- table(classifications$cluster,classifications$predicted)
  
  head(classifications_table)
  
  classifications_table <- as.matrix(classifications_table)
  
  M <- matrix(
    classifications_table, 
    nrow = length(unique(classifications$cluster)), 
    ncol = length(unique(classifications$predicted)), 
    dimnames = dimnames(classifications_table)
  )
  
  M[is.na(M)] <- 0
  
  #M[1:5,1:5]
  
  #
  relativise <- function(x) {
    return( (x - min(x)) / (max(x) - min(x)) )
  }
  
  
  M <- t(apply(M,1,relativise))
  
  #M[1:5,1:5]
  Heatmap(
    t(M),
    name = "% cells",
    cluster_rows = FALSE,
    col = viridis(10),
    clustering_method_columns = "complete",
    column_names_side = "top"
  )
  
  #OR
  
  M2 <- M[,apply(M,2,function(x){any(x>0.8)})] # 
  
  Heatmap(
    t(M2),
    name = "% cells",
    cluster_rows = FALSE,
    col = viridis(10),
    clustering_method_columns = "complete",
    column_names_side = "top"
  )
  
  M3 <- M2
  
  M3[M3 < 0.75] <- 0
  
  g <- graph_from_incidence_matrix(
    incidence = M3,
    directed = TRUE,
    mode = "in",
    weighted = TRUE,
    add.names = NULL
  )
  
  a=plot(
    g,
    vertex.color = cl_colors,
    vertex.size = 4,
    edge.size = 1,
    edge.arrow.size = 0.5,
    vertex.label.size = 0.1,
    layout = layout_as_bipartite(g)[,c(2,1)]
  )
  
  
  ht=Heatmap(
    unique(M3),
    name = "% cells",
    cluster_rows = FALSE,
    col = viridis(10),
    clustering_method_columns = "complete",
    column_names_side = "top"
  )
  
  tiff(paste0("./figures/cell_type/sizes_transferlabels_graph_",rc,".tiff"), 
       height = 40,width = 60,units = "cm",bg = NA, res=300, compression = "lzw")
  plot(ht)
  dev.off()
  
  
  pdf(paste0("./figures/cell_type/sizes_transferlabels_graph_",rc,".pdf"), wi = 50, he = 50)
  set.seed(123)
  plot(
    g,
    vertex.color = cl_colors,
    vertex.size = 4,
    edge.size = 1,
    edge.arrow.size = 0.5,
    vertex.label.size = 0.1,
    layout = layout_as_bipartite(g)[,c(2,1)]
  )
  dev.off()
  
  #anno_table
  write.csv(t(M3),paste0("./outputs/cell_type/cell_vote_",rc,".csv"))
  #by celltype
  tb=t(M3)
  ct_list=list()
  for (i in 1:nrow(t(M3))) {
    rm=sub(".*\\.", "", rownames(t(M3))[i])
    clt=colnames(tb)[which(tb[i,]== 1)]
    ct_list[[rm]]=c(ct_list[[rm]],clt)
  }
  
  #by cluster
  
  tb=M3
  ct_list2=list()
  for (i in 1:nrow(M3)) {
    rm=rownames(M3)[i]
    clt=colnames(tb)[which(tb[i,]== 1)]
    cm=sub(".*\\.", "", clt)
    ct_list2[[rm]]=unique(c(ct_list2$rm,cm))
    
  }
  ct_name=paste0("ct_list_",rc)
  assign(ct_name,ct_list2)
}

#####

#summary
#cell-type
for (rc in rc_list) {
  tname=paste0("ct_list_",rc)
  #for ref_cell
  df=get(tname)
  df=lapply(df, function(x) substr(x, 1, 5))
  
  # Convert to data frame
  df_out <- data.frame(
    cluster = names(df),
    label = sapply(df, function(x) {
      if (length(x) <= 2) {
        paste(x, collapse = ".")
      } else if (length(x) == 3|4) {
        paste(x[1:2], collapse = ".")
      } else {
        "unD"
      }
    }),
    stringsAsFactors = FALSE
  )
  
  df_out$label_numbered <- 
    paste0(df_out$cluster,"_",df_out$label)
  
  
  # Result
  print(df_out)
  
  snames=paste0("uniqu",rc)
  assign(snames,df_out)
  write.csv(df_out,paste0("./outputs/cell_type/enrichr_anno_",rc,".csv"))
  
  enrichr_anno=as.character(ss$merged_sub)
  enrichr_anno_type=as.character(ss$merged_sub)
  
  for (i in df_out$cluster) {
    enrichr_anno[enrichr_anno == i] = df_out$label[which(df_out$cluster == i)]
    enrichr_anno_type[enrichr_anno_type == i] = df_out$label_numbered[which(df_out$cluster == i)]
    print(i)
  }
  
  if (rc == "ref_cell1") {
    ss[["enrichr_anno_ref1"]]= enrichr_anno
    ss[["enrichr_anno_type_ref1"]]= enrichr_anno_type
    print(unique(enrichr_anno_type))
  }else{
    ss[["enrichr_anno_ref2"]]=enrichr_anno
    ss[["enrichr_anno_type_ref2"]]=enrichr_anno_type
    print(unique(enrichr_anno_type))
  }
}       

#
#sctype
source("./codes/scType.R")

sctype_scores=unique(sctype_scores)
write.csv(sctype_scores,paste0("./outputs/cell_type/sctype_scores_var",r.variable,".csv"))


#single-tag
# Function to convert phrase to abbreviation
abbreviate_label <- function(x) {
  words <- unlist(strsplit(x, " "))
  if (words[1] == "Unknown") {
    return(paste0("UnD"))
  } else if (words[1] == "Neuroblasts") {
    return(toupper(substring(words[1], 1, 6)))
  } else if (length(words) == 1) {
    return(toupper(substring(words[1], 1, 5)))
  } else if (length(words) == 2 ) {
    return(paste0(toupper(substring(words[1], 1, 4)), tolower(substring(words[2], 1, 1))))
  }else {
    return(paste0(toupper(substring(words[1], 1, 1)), 
                  toupper(substring(words[2], 1, 1)),
                  toupper(substring(words[3], 1, 1))
    ))
  }
}

# Apply function
abbreviated <- sapply(sctype_scores$type, abbreviate_label)

# Output result

sctype_scores$abb=abbreviated
sctype_scores$label=paste0(sctype_scores$cluster,"_",sctype_scores$abb)


df= sctype_scores %>% dplyr::group_by(cluster) %>% dplyr::summarise(abb = paste(unique(abb), collapse = "/"),
                                                                    label = paste0(cluster[1], "_", paste(unique(abb), collapse = "/")),
                                                                    .groups = "drop")

enrichr_anno=as.character(ss$merged_sub)
enrichr_anno_type=as.character(ss$merged_sub)

for (i in df$cluster) {
  enrichr_anno[enrichr_anno == i] = df$abb[which(df$cluster == i)]
  enrichr_anno_type[enrichr_anno_type == i] = df$label[which(df$cluster == i)]
  print(i)
}
ss[["enrichr_anno_ref3"]]= enrichr_anno
ss[["enrichr_anno_type_ref3"]]= enrichr_anno_type
print(unique(enrichr_anno_type))


top10.merged_mk= all.merged.markers %>%
  dplyr::filter(pct.1>0.3) %>%
  dplyr::filter(avg_log2FC>0.6) %>%
  dplyr::filter(p_val_adj <0.01) %>%
  group_by(cluster) %>%
  top_n(10, wt=(pct.1/pct.2))

type_table=data.frame(
                      "cluster"=unique(ss$merged_sub),
                      "ref1"=unique(ss$enrichr_anno_type_ref1),
                      "ref2"=unique(ss$enrichr_anno_type_ref2),
                      "ref3"=unique(ss$enrichr_anno_type_ref3))

top10mk=top10.merged_mk %>% dplyr::group_by(cluster) %>%
  dplyr::summarise(cluster = unique(cluster), 
                   topmk = paste(unique(gene), collapse = "/"),
                   .groups = "drop")


type_table_m = left_join(type_table,top10mk, by = "cluster" )
np.top100_mk=dplyr::filter(top10.merged_mk,gene %in% neuropep)
np.top100_mk
write.csv(type_table_m,paste0("./outputs/cell_type/cell_type_table_m.csv"))

dir.create(paste0("./figures/cell_type/mk_feature_plots_pre"))

for (i in 1:nrow(type_table_m)) {
  cl=type_table_m$cluster[i]
  mk=str_split(type_table_m$topmk[i],"\\/")[[1]]
  if (length(mk)>1) {
    print(i)
    print(cl)
    tiff(paste0("./figures/cell_type/mk_feature_plots_pre/",cl,"_",r.variable,"_mk_FeaturePlot.tiff"),
         width = (10*round(length(mk)/2))+10,
         height = (10*round(length(mk)/(round(length(mk)/2)))),
         units = "cm", res = 300,compression = "lzw",bg = NA)
    p=FeaturePlot(ss,
                  features = mk,
                  pt.size = 0.7, cols = c("grey80","red"),
                  reduction=paste0(umap),
                  order = T,min.cutoff = 0.8,
                  ncol=round(length(mk)/2)) & NoAxes()& NoLegend()
    p0=CellDimPlot(ss, seed = 0,group.by = "merged_sub",
                  reduction=umap, 
                  cells.highlight = names(ss$merged_sub[which(ss$merged_sub == cl)])) 
    
    print((p|p0)+plot_layout(widths = c(2, 1)))
    dev.off()
  }else if(length(mk)==1){
    print("no marker")}
}


#visualization
tiff(paste0("./figures/cell_type/",r.variable,"enrichr_annotations.tiff"),
     width = 60,height = 60,units = "cm", res = 300,compression = "lzw", bg = NA)

p1=CellDimPlot(
  srt = ss, group.by = "enrichr_anno_type_ref1", seed = 0,
  reduction = umap, theme_use = "theme_blank",legend.position = "none",
  label = T,label_insitu = T,label_repel = T
)
p2=CellDimPlot(
  srt = ss, group.by = "enrichr_anno_type_ref2", seed = 0,
  reduction = umap, theme_use = "theme_blank",legend.position = "none",
  label = T,label_insitu = T,label_repel = T
)
p3=CellDimPlot(
  srt = ss, group.by = "enrichr_anno_type_ref3", seed = 0,
  reduction = umap, theme_use = "theme_blank",legend.position = "none",
  label = T,label_insitu = T,label_repel = T
)

p4=CellDimPlot(
  srt = ss, group.by = "enrichr_anno_ref1", seed = 0,
  reduction = umap, theme_use = "theme_blank",legend.position = "none",
  label = T,label_insitu = T,label_repel = T
)

p5=CellDimPlot(
  srt = ss, group.by = "enrichr_anno_ref2", seed = 0,
  reduction = umap, theme_use = "theme_blank",legend.position = "none",
  label = T,label_insitu = T,label_repel = T
)
p6=CellDimPlot(
  srt = ss, group.by = "enrichr_anno_ref3", seed = 0,
  reduction = umap, theme_use = "theme_blank",legend.position = "none",
  label = T,label_insitu = T,label_repel = T
)

print((p1|p2|p3)/(p4|p5|p6))

dev.off()

tiff(paste0("./figures/cell_type/",r.variable,"enrichr_annotations_label_free.tiff"),
     width = 60,height = 20,units = "cm", res = 300,compression = "lzw", bg = NA)
p1=CellDimPlot(
  srt = ss, group.by = "enrichr_anno_type_ref1", seed = 0,
  reduction = umap, theme_use = "theme_blank",
  label = F,label_insitu = T,label_repel = T
)
p2=CellDimPlot(
  srt = ss, group.by = "enrichr_anno_type_ref2", seed = 0,
  reduction = umap, theme_use = "theme_blank",
  label = F,label_insitu = T,label_repel = T
)
p3=CellDimPlot(
  srt = ss, group.by = "enrichr_anno_type_ref3", seed = 0,
  reduction = umap, theme_use = "theme_blank",
  label = F,label_insitu = T,label_repel = T
)

p4=CellDimPlot(
  srt = ss, group.by = "enrichr_anno_ref1", seed = 0,
  reduction = umap, theme_use = "theme_blank",
  label = F,label_insitu = T,label_repel = T
)

p5=CellDimPlot(
  srt = ss, group.by = "enrichr_anno_ref2", seed = 0,
  reduction = umap, theme_use = "theme_blank",
  label = F,label_insitu = T,label_repel = T
)
p6=CellDimPlot(
  srt = ss, group.by = "enrichr_anno_ref3", seed = 0,
  reduction = umap, theme_use = "theme_blank",
  label = F,label_insitu = T,label_repel = T
)

print((p1|p2|p3)/(p4|p5|p6))

dev.off()


tiff(paste0("./figures/cell_type/",r.variable,"merged_sube.tiff"),
     width = 60,height = 40,units = "cm", res = 300,compression = "lzw", bg = NA)
p1=CellDimPlot(
  srt = ss, group.by = "merged_sub", seed = 0,
  reduction = "umap", theme_use = "theme_blank",
  label = T,label_insitu = T,label_repel = T
)
p2=CellDimPlot(
  srt = ss, group.by = "merged_sub", seed = 0,
  reduction = "umap", theme_use = "theme_blank",
  label = F,label_insitu = T,label_repel = T
)

p3=CellDimPlot(
  srt = ss, group.by = "merged_sub", seed = 0,
  reduction = umap, theme_use = "theme_blank",
  label = T,label_insitu = T,label_repel = T
)
p4=CellDimPlot(
  srt = ss, group.by = "merged_sub", seed = 0,
  reduction = umap, theme_use = "theme_blank",
  label = F,label_insitu = T,label_repel = T
)
print((p1|p2)/(p3|p4))

dev.off()
#merge_celltype_annotation
#saveRDS(ss,paste0("./data/rds/step4_var",r.variable,".rds"))


##manual correction
type_table_m.v2= read.csv(paste0("./outputs/cell_type/cell_type_table_m_modified.csv"), header = T, row.names = 1)

merged_sub.anno= ss$merged_sub
merged_sub.anno_type= ss$merged_sub
merged_sub.anno_type_label= ss$merged_sub
merged_sub_numb= ss$merged_sub

for (i in type_table_m.v2$cluster) {
  merged_sub.anno[merged_sub.anno == i] = type_table_m.v2$broader[which(type_table_m.v2$cluster == i)]
  merged_sub.anno_type[merged_sub.anno_type == i] = type_table_m.v2$numbered[which(type_table_m.v2$cluster == i)]
  merged_sub.anno_type_label[merged_sub.anno_type_label == i] = type_table_m.v2$numbered_mk[which(type_table_m.v2$cluster == i)]
  merged_sub_numb[merged_sub_numb == i] = type_table_m.v2$cl_numb[which(type_table_m.v2$cluster == i)]
  print(i)
}

ss[["merged_sub.anno"]]=merged_sub.anno
ss[["merged_sub.anno_type"]]=merged_sub.anno_type
ss[["merged_sub.anno_type_label"]]=merged_sub.anno_type_label
ss[["merged_sub_numb"]]=merged_sub_numb

ss$merged_sub.anno_type_label=factor(ss$merged_sub.anno_type_label,levels = unique(type_table_m.v2$numbered_mk))
ss$merged_sub.anno_type=factor(ss$merged_sub.anno_type,levels = unique(type_table_m.v2$numbered))
ss$merged_sub.anno=factor(ss$merged_sub.anno,levels = unique(type_table_m.v2$broader))
ss$merged_sub_numb=factor(ss$merged_sub_numb,levels = unique(type_table_m.v2$cl_numb))


a=VlnPlot(ss,group.by = "merged_sub.anno_type",split.by = "orig.ident",
          cols = rev(c("blue","lightgreen","orange","red")),
          features = c("nCount_RNA","nFeature_RNA"),ncol = 1,
          pt.size = 0)
a


tiff(paste0("./figures/cell_distribution_by_cluster_merged_sub",r.variable,".tiff"),width = 80,height = 30,units = "cm", 
     res = 300,compression = "lzw")
print(a)
dev.off()

#### cyto color

#broad color
unique(type_table_m.v2$broader)
broader_col= c(
  "#ff8f00","darkgreen","#1F78B4","#ea80fc","purple","#cfd8dc","#f50057"
)
names(broader_col)=unique(type_table_m.v2$broader)

type_table_m.v2$broader_col=type_table_m.v2$broader
type_table_m.v2$broader_col[type_table_m.v2$broader_col== names(broader_col)[1]]=broader_col[1]
type_table_m.v2$broader_col[type_table_m.v2$broader_col== names(broader_col)[2]]=broader_col[2]
type_table_m.v2$broader_col[type_table_m.v2$broader_col== names(broader_col)[3]]=broader_col[3]
type_table_m.v2$broader_col[type_table_m.v2$broader_col== names(broader_col)[4]]=broader_col[4]
type_table_m.v2$broader_col[type_table_m.v2$broader_col== names(broader_col)[5]]=broader_col[5]
type_table_m.v2$broader_col[type_table_m.v2$broader_col== names(broader_col)[6]]=broader_col[6]
type_table_m.v2$broader_col[type_table_m.v2$broader_col== names(broader_col)[7]]=broader_col[7]

#cyt_color
table(type_table_m.v2$broader)

ref_colors <- broader_col

lighten_steps <- function(color, n) {
  lighten_seq <- seq(0, 0.6, length.out = n)  
  lighten(color, amount = lighten_seq)
}

type_table_m.v2$clcol=type_table_m.v2$broader
type_table_m.v2[which(type_table_m.v2$broader == names(broader_col)[1]),"clcol"]=lighten_steps(ref_colors[1],sum(type_table_m.v2$broader== names(broader_col)[1])) #"#"
type_table_m.v2[which(type_table_m.v2$broader == names(broader_col)[2]),"clcol"]=lighten_steps(ref_colors[2],sum(type_table_m.v2$broader== names(broader_col)[2])) #"#"
type_table_m.v2[which(type_table_m.v2$broader == names(broader_col)[3]),"clcol"]=lighten_steps(ref_colors[3],sum(type_table_m.v2$broader== names(broader_col)[3])) #"#"
type_table_m.v2[which(type_table_m.v2$broader == names(broader_col)[4]),"clcol"]=lighten_steps(ref_colors[4],sum(type_table_m.v2$broader== names(broader_col)[4])) #"#"
type_table_m.v2[which(type_table_m.v2$broader == names(broader_col)[5]),"clcol"]=lighten_steps(ref_colors[5],sum(type_table_m.v2$broader== names(broader_col)[5])) #"#"
type_table_m.v2[which(type_table_m.v2$broader == names(broader_col)[6]),"clcol"]=lighten_steps(ref_colors[6],sum(type_table_m.v2$broader== names(broader_col)[6])) #"#"
type_table_m.v2[which(type_table_m.v2$broader == names(broader_col)[7]),"clcol"]=lighten_steps(ref_colors[7],sum(type_table_m.v2$broader== names(broader_col)[7])) #"#"

write.csv(type_table_m.v2,paste0("./outputs/cell_type/cell_type_table_m_modified.csv"))

#marker UMAP
head(type_table_m.v2)
dir.create(paste0("./figures/cell_type/mk_feature_plots"))
for (i in 1:nrow(type_table_m.v2)) {
  cl=gsub("\\/",".",type_table_m.v2$numbered[i])
  mk=str_split(type_table_m.v2$topmk[i],"\\/")[[1]]
  if (length(mk)>1) {
    print(i)
    print(cl)
    tiff(paste0("./figures/cell_type/mk_feature_plots/",cl,"_",r.variable,"_mk_FeaturePlot.tiff"),
         width = (10*round(length(mk)/2)),
         height = 10*round(length(mk)/(round(length(mk)/2))),
         units = "cm", res = 300,compression = "lzw",bg = NA)
    p=FeaturePlot(ss,
                  features = mk,
                  pt.size = 0.7, cols = c("grey80","red"),
                  reduction=paste0(umap),
                  order = T,min.cutoff = 0.8,
                  ncol=round(length(mk)/2)) & NoAxes()& NoLegend()
    
    print(p)
    dev.off()
  }else if(length(mk)==1){
    print("no marker")}
}

#RNA map
for (i in 1:nrow(type_table_m.v2)) {
  cl=gsub("\\/",".",type_table_m.v2$numbered[i])
  mk=str_split(type_table_m.v2$topmk[i],"\\/")[[1]]
  if (length(mk)>1) {
    print(i)
    print(cl)
    tiff(paste0("./figures/cell_type/mk_feature_plots/RNA_only_",cl,"_",r.variable,"_mk_FeaturePlot.tiff"),
         width = (10*round(length(mk)/2)),
         height = 10*round(length(mk)/(round(length(mk)/2))),
         units = "cm", res = 300,compression = "lzw",bg = NA)
    p=FeaturePlot(ss,
                  features = mk,
                  pt.size = 0.7, cols = c("grey80","red"),
                  reduction="umap",
                  order = T,min.cutoff = 0.8,
                  ncol=round(length(mk)/2)) & NoAxes()& NoLegend()
    
    print(p)
    dev.off()
  }else if(length(mk)==1){
    print("no marker")}
}
#### cluster color setting

ss4000_ctypes= ss[["merged_sub.anno_type"]]
colnames(ss4000_ctypes)="ctype"
ss4000_ctypes$simp_ctype=factor(ss$merged_sub.anno,levels = sort(unique(type_table_m.v2$broader)))

####cluster meta
ss4000_ctypes
Idents(ss) <- ss$merged_sub.anno_type
#meta col
#1 ctype_col

ss4000_ctypes$ctype_col=translate_ids(ss4000_ctypes$ctype, type_table_m.v2[,c("numbered","clcol")])

#2simp_ctype_col

ss4000_ctypes$simp_ctype_col=translate_ids(ss4000_ctypes$simp_ctype, type_table_m.v2[,c("broader","broader_col")])

ss[["ctype_col"]]=ss4000_ctypes$ctype_col
ss[["simp_ctype_col"]]=ss4000_ctypes$simp_ctype_col



###dpt
library(destiny)
dm <- DiffusionMap(Embeddings(ss, "harmony")[,1:149])
dpt <- DPT(dm)
ss$dpt <- rank(dpt$dpt)

FeaturePlot(
  ss, features = "dpt", split.by = "orig.ident",
  reduction = umap
)& NoLegend()

#hierarchical clustering 
#simple test
#ss <- BuildClusterTree(object = ss,reduction = "harmony")
#a=PlotClusterTree(object = ss)

saveRDS(ss,paste0("./data/rds/norm_step4_withUnD_var",r.variable,".rds"))

source("./codes/04_c_hclust_cluster_06_2025.R")


#subsetting_removing UnD clusters
'%!in%' <- function(x,y)!('%in%'(x,y))
unD= unique(ss$merged_sub.anno_type)[grep("UnD",unique(ss$merged_sub.anno_type) )]
ss1=subset(ss, merged_sub.anno_type %!in%  unD )
type_table_m.v3= type_table_m.v2[- c(grep(pattern = "UnD", type_table_m.v2$numbered)),]
ss1$merged_sub.anno_type=factor(ss1$merged_sub.anno_type, levels = type_table_m.v3$numbered)

ss4000_ctypes=ss4000_ctypes[- c(grep(pattern = "UnD", ss4000_ctypes$ctype)),]

write.csv(
  ss4000_ctypes,
  file = paste0("./outputs/cell_type/ss4000_ctypes.csv"),
  row.names = T,
  quote = FALSE
)

simcol=type_table_m.v2$broader_col
names(simcol)=type_table_m.v2$broader
p0=DimPlot(ss,group.by = "merged_sub.anno", seed = 0,
           reduction = umap,cols = simcol,alpha = 0.8,
           label = F, label.size =3,stroke.size = 0.5,
           order=T
)+ ggtitle(NULL) & NoLegend()

cycol=type_table_m.v2$clcol
names(cycol)=type_table_m.v2$numbered
p01=DimPlot(ss,group.by = "merged_sub.anno_type", seed = 0,
            reduction = umap,cols = cycol,alpha = 0.8,
            label = F, label.size =3,stroke.size = 0.5,
            order=T
)+ ggtitle(NULL) & NoLegend()
  

p1=DimPlot(ss1,group.by = "merged_sub.anno", seed = 0,
           reduction = umap,cols = simcol,alpha = 0.8,
           label = F, label.size =3,stroke.size = 0.5,
           order=T
)+ ggtitle(NULL) & NoLegend()
  

p2=DimPlot(ss1,group.by = "merged_sub.anno_type", seed = 0,
           reduction = umap,cols = cycol,alpha = 0.8,
           label = F, label.size =3,stroke.size = 0.5,
           order=T
)+ ggtitle(NULL) & NoLegend()


p3=DimPlot(ss1,group.by = "merged_sub.anno", seed = 0,
           reduction = umap,cols = simcol,alpha = 0.8,
           label = T, label.size =3,stroke.size = 0.5,
           order=T
)+ ggtitle(NULL) & NoLegend()
  

p4=DimPlot(ss1,group.by = "merged_sub.anno_type", seed = 0,
           reduction = umap,cols = cycol,alpha = 0.8,
           label = T, label.size =3,stroke.size = 0.5,
           order=T
)+ ggtitle(NULL) & NoLegend()

p5= DimPlot(ss1,group.by = "seurat_clusters", seed = 0,
            reduction = umap,alpha = 0.8,
            label = T, label.size =3,stroke.size = 0.5,
            order=T
)+ ggtitle(NULL) & NoLegend()
  
  
p6=DimPlot(ss1,group.by = "merged_sub", seed = 0,
           reduction = umap,cols = cycol,alpha = 0.8,
           label = T, label.size =3,stroke.size = 0.5,
           order=T
)+ ggtitle(NULL) & NoLegend()
  

pdf(paste0("./figures/cell_type/",r.variable,"submerge_annotation.pdf"),
     width = 40,height = 40)
print((p1|p2)/(p3|p4))

dev.off()

pdf(paste0("./figures/cell_type/",r.variable,"submerge_annotation_UnD.pdf"),
     width = 40,height = 20)
print((p0|p01))

dev.off()

pdf(paste0("./figures/cell_type/",r.variable,"submerg.pdf"),
     width = 40,height = 20)
print(p5|p6)

dev.off()


p7=DimPlot(ss1,group.by = "merged_sub.anno", seed = 0,
           reduction = umap,cols = simcol,alpha = 0.8,
           label = F, label.size =3,stroke.size = 0.5
)+ ggtitle(NULL)+theme(text = element_text(size = 12),axis.title =element_text(size = 10) )

numcol=type_table_m.v3$clcol
names(numcol)=type_table_m.v3$cl_numb
library(shadowtext)
p8=DimPlot(ss1,group.by = "merged_sub_numb", seed = 0,
           reduction = umap,cols = numcol,alpha = 0.5,
           label = T,label.size=3,label.color = c("black"), stroke.size = 0.5,
           order=T
           )+ ggtitle(NULL) +theme(text = element_text(size = 12),axis.title =element_text(size = 10) )& NoLegend()

p9=DimPlot(ss1,group.by = "merged_sub_numb", seed = 0,
           reduction = "umap",cols = numcol,alpha = 0.5,
           label = T,label.size=3,label.color = c("black"), stroke.size = 0.5,
           order=T
)+ ggtitle(NULL)+theme(text = element_text(size = 12),axis.title =element_text(size = 10) )& NoLegend()  


pdf(paste0("./figures/cell_type/npcs_umap",r.variable,".pdf"),
     width = 18,height = 5)
print(p7|p8|p9)
dev.off()

pdf(paste0("./figures/cell_type/npcs_umap_sp",r.variable,".pdf"),
    width = 5,height = 18)
p10=DimPlot(ss1,group.by = "merged_sub_numb", seed = 0,split.by = "orig.ident",
        reduction = umap,cols = numcol,
        label = F,label.size=3,label.color = c("black"), stroke.size = 0.5,ncol = 1,
        order=T
)+ ggtitle(NULL) +theme(text = element_text(size = 12),axis.title =element_text(size = 10) )& NoLegend()
print(p10)
dev.off()
## cluster marker summary

all.merge.anno.markers= FindAllMarkers(
  ss1,assay = "RNA",
  group.by = "merged_sub.anno_type",
  logfc.threshold = 0.1,
  test.use = "wilcox",#"DESeq2","MAST","wilcox"
  #slot = "data",
  min.pct = 0.01,
  min.diff.pct = -Inf,
  node = NULL,
  verbose = TRUE,
  only.pos = T,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  mean.fxn = NULL,
  fc.name = NULL,
  base = 2,
  return.thresh = 0.01,
  densify = FALSE
)
#top100s
top100_sub_mk= all.merge.anno.markers %>% 
  dplyr::filter(p_val_adj <0.05) %>%
  group_by(cluster) %>%
  top_n(100, wt=avg_log2FC)

write.csv(all.merge.anno.markers,paste0("./outputs/cell_type/all.sub.markers_",
                                        "var",r.variable,
                                        "_merged_sub.anno.csv"))

top5_sub_mk= all.merge.anno.markers %>% 
  dplyr::filter(avg_log2FC>1) %>%
  dplyr::filter(p_val_adj <0.05) %>%
  group_by(cluster) %>%
  top_n(5, wt=avg_log2FC)

#all_sig
all_sig_mk= all.merge.anno.markers %>% 
  dplyr::filter(avg_log2FC>1) %>%
  dplyr::filter(p_val_adj <0.05) %>%
  dplyr::filter(pct.1 >0.5) %>%
  group_by(cluster) 
ht1 <- GroupHeatmap(ss1,assay = "RNA",
                    features = all_sig_mk$gene, feature_split = all_sig_mk$cluster, 
                    row_title = "",heatmap_palcolor = c("white","#fffef5","#fffde7","#fff9ab","#fff5d3","#fd8d3c"),
                    group.by = "merged_sub.anno_type",#split.by = "orig.ident", #group.by = "sub_d25_res0.5"
                    flip = T,column_title_rot = 45,group_palette = "simspec",
                    #cell_annotation = c(nps),#,"Phase","foxp2","foxp4", "oxt"), 
                    #cell_annotation_palette = c(rep("Paired",length(nps))),#, Dark2"Paired", "Paired", "Paired"),
                    #cell_annotation_params = list(height = unit(10, "mm")),
                    #feature_annotation = c("TF","dGCprimed"),label_size = 2,
                    #feature_annotation_palcolor = list(c("gold", "steelblue"),c("red")),
                    add_dot = F, add_bg = T,border = T, nlabel = 0, show_row_names = TRUE
)


pdf(paste0("./figures/cell_type/all_orig_mk_r4000_mergedsub.pdf"),
     width = 90,height = 20)
print(ht1$plot)
dev.off()

#for top5
ht1 <- GroupHeatmap(ss1,assay = "RNA",
                    features = top5_sub_mk$gene, feature_split = top5_sub_mk$cluster, 
                    row_title = "",
                    group.by = "merged_sub.anno_type",#split.by = "orig.ident", #group.by = "sub_d25_res0.5"
                    heatmap_palette = "YlGn",flip = T,column_title_rot = 30,group_palette = "simspec",
                    #cell_annotation = c(nps),#,"Phase","foxp2","foxp4", "oxt"), 
                    #cell_annotation_palette = c(rep("Paired",length(nps))),#, Dark2"Paired", "Paired", "Paired"),
                    #cell_annotation_params = list(height = unit(10, "mm")),
                    #feature_annotation = c("TF","dGCprimed"),label_size = 2,
                    #feature_annotation_palcolor = list(c("gold", "steelblue"),c("red")),
                    add_dot = F, add_bg = T,border = T, nlabel = 0, show_row_names = TRUE
)


pdf(paste0("./figures/cell_type/top5_orig_mk_r3500_mergedsub.pdf"),
     width = 90,height = 20)
print(ht1$plot)
dev.off()



# cell type data save -----------------------------------------------------
'%!in%'=function(x,y)!('%in%'(x,y))
#ss1$merged_sub.anno_type=factor(ss1$merged_sub.anno_type, levels =type_table_m.v3$numbered)

ss_col <- 
  setNames(
    ss1$ctype_col,
    ss1$merged_sub.anno_type
  )

#r consensus_UMAP,  fig.width = 14, fig.height = 6, warning = FALSE, message = FALSE}

pdf(paste0("./figures/cell_type/",r.variable,"nps_sub.pdf"),
     width = 25,height = 20)

a=DimPlot(
  ss1,
  reduction = umap,
  group.by = "merged_sub.anno_type",
  cols = ss_col,alpha = 0.5,
  pt.size = 0.7,
  raster = FALSE
) + NoAxes()

print(a)
dev.off()

#saveRDS(ss1,paste0("./data/rds/norm_step4_var",r.variable,".rds"))

save(
  ss4000_ctypes,
  taxa_order,
  all.merge.anno.markers,
  top100_sub_mk,
  file = paste0("./data/rda/ss1_step4_before_scores.rda")
)

###
load( file = paste0("./data/rda/ss1_step4_before_scores.rda"))
#milo, test differential abundance 
source("./codes/04_b_MiloR_run.R")
#####
##Identification of LD-responsive cell clusters


#1.AUC and Ucell score calculation
#activity markers, early responding genes: egr+ cells (fosaa, fosab, ...)
#GC-receptors: nr3c1,nr3c2
#GC-response marker gene: fkbp5

source("./codes/04_a_AUC_UCell_score_calculation.R")


######## 

# gene Activity ----------------------------------------------------------------

#This is to compare with what we obtained from the scATAC-seq geneActivity. The intersection of the two will be used to select genomic regions of interest.
## scRNA-seq broad type markers

#This is to compare with what we obtained from the scATAC-seq geneActivity. The intersection of the two will be used to select genomic regions of interest.

Idents(ss1) <- ss1$merged_sub.anno_type
DefaultAssay(ss1)="ATAC"

vio_scatac_feat_peaks <- 
  VlnPlot(ss1, assay = "ATAC", 
          features = c("nFeature_ATAC"),
          group.by = "merged_sub.anno_type",cols = unique(ss1$ctype_col), pt.size = 0)

set.seed(456)
feat_subspl <- vio_scatac_feat_peaks$data %>%
  group_by(ident) %>%
  sample_n(250, replace = TRUE) %>%
  ungroup()

library(colorspace)

vio_scatac_feat_peaks <- vio_scatac_feat_peaks + geom_jitter(
  mapping = aes(color = ident), data = feat_subspl,
  position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9))+
  scale_color_manual(breaks = ss1$merged_sub.anno_type, 
                     values = darken(alpha(ss1$ctype_col, .3),.5)
  )+
  theme(legend.position = "none")


vio_scatac_umi_peaks <- 
  VlnPlot(ss1, features = c("nFeature_ATAC"),
          group.by = "merged_sub.anno_type",cols = ss1$ctype_col, pt.size = 0)

set.seed(456)
umi_subspl <- vio_scatac_umi_peaks$data %>%
  group_by(ident) %>%
  sample_n(250, replace = TRUE) %>%
  ungroup()

vio_scatac_umi_peaks <- vio_scatac_umi_peaks + geom_jitter(
  mapping = aes(color = ident), data = umi_subspl,
  position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9))+
  scale_color_manual(breaks = ss1$merged_sub.anno_type, values = darken(alpha(ss1$ctype_col, .3),.5))+
  theme(legend.position = "none")


pdf(paste0("./figures/cell_type/supp_1_1_CD.pdf"), width = 10, height = 10)
cowplot::plot_grid(
  vio_scatac_feat_peaks,vio_scatac_umi_peaks, labels = c("A","B"),ncol = 1
)
dev.off()

vio_scatac_feat_genact <- 
  VlnPlot(ss1, features = c("nFeature_RNA"),
          group.by = "merged_sub.anno_type",cols = ss1$ctype_col, pt.size = 0)

set.seed(456)
feat_subspl <- vio_scatac_feat_genact$data %>%
  group_by(ident) %>%
  sample_n(250, replace = TRUE) %>%
  ungroup()

vio_scatac_feat_genact <- vio_scatac_feat_genact + geom_jitter(
  mapping = aes(color = ident), data = feat_subspl,
  position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9))+
  scale_color_manual(breaks = ss1$merged_sub.anno_type, values = darken(alpha(ss1$ctype_col, .3),.5))+
  theme(legend.position = "none")

#RNA
vio_scatac_umi_genact <- 
  VlnPlot(ss1, features = c("nCount_RNA"),
          group.by = "merged_sub.anno_type",cols = ss1$ctype_col, pt.size = 0)

set.seed(456)
umi_subspl <- vio_scatac_umi_genact$data %>%
  group_by(ident) %>%
  sample_n(250, replace = TRUE) %>%
  ungroup()

vio_scatac_umi_genact <- vio_scatac_umi_genact + geom_jitter(
  mapping = aes(color = ident), data = umi_subspl,
  position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9))+
  scale_color_manual(breaks = ss1$merged_sub.anno_type, values = darken(alpha(ss1$ctype_col, .3),.5))+
  theme(legend.position = "none")


pdf(paste0("./figures/cell_type/supp_1_1_EF.pdf"), width = 10, height = 10)
cowplot::plot_grid(
  vio_scatac_feat_genact,vio_scatac_umi_genact, labels = c("A","B"),ncol = 1
)
dev.off()


# gene score --------------------------------------------------------------

####gene score
DefaultAssay(ss1)="RNA"

##gene score
mat = ss1@assays$RNA$data
all_genes = rownames(mat)

scores <- list()

for(i in str_sort(unique(top100_sub_mk$cluster),numeric = T)){
  marks = top100_sub_mk$gene[top100_sub_mk$cluster == i]
  marks = marks[marks %in% rownames(mat)]
  scores[[i]] <- 
    gene_score(
      x = mat, gene_set = marks,
      gene_pool = all_genes, remove_set_from_pool = TRUE,
      fraction = 0.05
    )
  scores[[i]] <- relativise(scores[[i]])
  message(i)
}

df_scores <-
  cbind(
    cell = colnames(ss1), # cells
    clu = ss1$merged_sub.anno_type, # cluster
    as.data.frame(scores)
  )

list_plots <- list()

for(i in names(scores)){
  j = gsub("/",".",i)
  j = gsub("\\+",".",j)
  message(i)
  list_plots[[i]] <-
    ggplot(
      data = df_scores,
      mapping =
        aes(x = clu,
            y = .data[[paste0("X",j)]],
            fill = clu,
            color = clu)) +
    geom_boxplot(outlier.colour = rgb(0.1,0.1,0.1,0.1), outlier.size = .4)+
    theme_classic()+
    scale_fill_manual(values = ss1$ctype_col)+
    scale_color_manual(values = darken(ss1$ctype_col,0.5))+
    ggtitle(paste0(i," RNA markers geneAct score"))+
    theme(axis.text.x = element_text(angle = 90,size = 6), legend.position = "none", plot.title = element_text(size = 8, face = "bold"))+
    coord_cartesian(ylim=c(0,1))
}

pdf(paste0("./figures/cell_type/markers_broad_geneAct_score.pdf"),width=15, height = 50)
cowplot::plot_grid(plotlist = list_plots, labels = LETTERS[1:9], ncol = 3)
dev.off()


#save as an object:

saveRDS(
  ss1,
  file = paste0("./data/rds/step4_var",r.variable,".rds")
)


## Plotting the number of cells per cluster
#modified code from from "Alberto Perez-Posada @apposada"
# create a data frame with the number of cells per cluster per library:

df <- 
  data.frame(
    table(
      ss1$orig.ident, factor(ss1$merged_sub.anno_type,levels = str_sort(unique(ss1$merged_sub.anno_type),numeric = T) )
    )
  )

# Rename the columns of the data frame
colnames(df) <- 
  c("Library", "Cluster", "Count")

# Normalize the counts in each cluster
df$logCount <- 
  log(df$Count,10) / log(tapply(df$Count, df$Cluster, sum)[df$Cluster],10)

df$CountNorm <- df$Count / tapply(df$Count, df$Cluster, sum)[df$Cluster]

# Create a stacked barplot raw number
ncells_bar <- ggplot(df, aes(x = Cluster, y = Count, fill = Library)) +
  geom_col(position = "stack") +
  scale_fill_brewer(
    type = "div",
    palette = "Spectral",
    direction = 1,
    aesthetics = "fill"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Number of cells") +
  guides(fill = FALSE)

# Create a stacked barplot log number
logncells_bar <- ggplot(df, aes(x = Cluster, y = logCount, fill = Library)) +
  geom_col(position = "stack") +
  scale_fill_brewer(
    type = "div",
    palette = "Spectral",
    direction = 1,
    aesthetics = "fill"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("log(Number of cells)") +
  guides(fill = FALSE)

# Create a stacked barplot norm number
normncells_bar <- ggplot(df, aes(x = Cluster, y = CountNorm, fill = Library)) +
  geom_col(position = "stack") +
  scale_fill_brewer(
    type = "div",
    palette = "Spectral",
    direction = 1,
    aesthetics = "fill"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Proportion of cells")

# Create a grid of the aligned plots
grid <- plot_grid(
  ncells_bar, logncells_bar, normncells_bar,
  nrow = 3,
  align = "v",
  axis = "tb",
  labels = c("A", "B", "C")
)

# Plot ncells per cluster
grid

tiff(paste0("./figures/cell_type/var",r.variable,"sub_merged__cell_numbers_ratio.tiff"),
     width = 30,height = 30,units = "cm", res = 300,compression = "lzw",bg = NA)
print(grid)

dev.off()

###

df <- 
  data.frame(
    table(
      ss1$orig.ident, factor(ss1$merged_sub.anno_type,levels = str_sort(unique(ss1$merged_sub.anno_type),numeric = T) )
    )
  )


# Rename the columns of the data frame
colnames(df) <- 
  c("Library", "Cluster", "Count")

# Normalize the counts in each cluster
df$logCount <- 
  log(df$Count,10) / log(tapply(df$Count, df$Cluster, sum)[df$Cluster],10)

df$CountNorm <- df$Count / tapply(df$Count, df$Cluster, sum)[df$Cluster]

df$Cluster=factor(df$Cluster, levels = type_table_m.v3$numbered)

 
# Create a stacked barplot raw number
ncells_bar <- ggplot(df, aes(x = Cluster, y = Count, fill = Library)) +
  geom_col(position = "stack") +
  scale_fill_brewer(
    type = "div",
    palette = "Spectral",
    direction = 1,
    aesthetics = "fill"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank(),strip.background = element_rect(colour="gray")) +
  ylab("Number of cells") +
  guides(fill = FALSE)

# Create a stacked barplot log number
logncells_bar <- ggplot(df, aes(x = Cluster, y = logCount, fill = Library)) +
  geom_col(position = "stack") +
  scale_fill_brewer(
    type = "div",
    palette = "Spectral",
    direction = 1,
    aesthetics = "fill"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank(),strip.background = element_rect(colour="gray")) +
  ylab("log(Number of cells)") +
  guides(fill = FALSE)

# Create a stacked barplot norm number
normncells_bar <- ggplot(df, aes(x = Cluster, y = CountNorm, fill = Library)) +
  geom_col(position = "stack") +
  scale_fill_brewer(
    type = "div",
    palette = "Spectral",
    direction = 1,
    aesthetics = "fill"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank(),strip.background = element_rect(colour="gray")) +
  ylab("Proportion of cells")

# Create a grid of the aligned plots
grid <- plot_grid(
  ncells_bar, logncells_bar, normncells_bar,
  nrow = 3,
  align = "v",
  axis = "tb"#,
  #labels = c("A", "B", "C")
)

# Plot ncells per cluster


tiff(paste0("./figures/cell_type/var",r.variable,"sub_merged__cell_numbers_ratio_sorted.tiff"),
     width = 30,height = 30,units = "cm", res = 300,compression = "lzw",bg = NA)
print(grid)

dev.off()


#####
#for the fig1


df <- 
  data.frame(
    table(
      ss$orig.ident, factor(ss$merged_sub.anno_type)
    )
  )

# Rename the columns of the data frame
colnames(df) <- 
  c("Library", "Cluster", "Count")

df$cl.n=sapply(str_split(df$Cluster,"_"),function (x){x[[1]]})

# Normalize the counts in each cluster
df$logCount <- 
  log(df$Count,10) / log(tapply(df$Count, df$Cluster, sum)[df$Cluster],10)

df$CountNorm <- df$Count / tapply(df$Count, df$Cluster, sum)[df$Cluster]
df$Cluster=gsub("\\/",".",df$Cluster)
df$Cluster=factor(df$Cluster, levels = type_table_m.v2$numbered)

##reorder colonm by hclust result
taxa_order_tb=data.frame("label"=taxa_order) 

tmp=join_by("label"=="numbered_mk")

taxa_order_tb=unique(left_join(taxa_order_tb,type_table_m.v2,by=tmp))

df$Cluster=factor(df$Cluster, levels =taxa_order_tb$numbered )

# Create a stacked barplot raw number
ncells_bar <- ggplot(df, aes(x = Cluster, y = Count, fill = Library)) +
  geom_col(position = "stack") +
  scale_fill_brewer(
    type = "div",
    palette = "Spectral",
    direction = 11,
    aesthetics = "fill"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.2),
        panel.background = element_blank(),strip.background = element_rect(colour="gray")) +
  ylab("Number of cells") +
  guides(fill = FALSE)

# Create a stacked barplot log number
logncells_bar <- ggplot(df, aes(x = Cluster, y = logCount, fill = Library)) +
  geom_col(position = "stack") +
  scale_fill_brewer(
    type = "div",
    palette = "Spectral",
    direction = 11,
    aesthetics = "fill"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.2),
        panel.background = element_blank(),strip.background = element_rect(colour="gray")) +
  ylab("log(Number of cells)") +
  guides(fill = FALSE)

# Create a stacked barplot norm number
normncells_bar <- ggplot(df, aes(x = Cluster, y = CountNorm, fill = Library)) +
  geom_col(position = "stack") +
  scale_fill_brewer(
    type = "div",
    palette = "Spectral",
    direction = 1,
    aesthetics = "fill"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.2),
        panel.background = element_blank(),strip.background = element_rect(colour="gray")) +
  ylab("Proportion of cells")

# Create a grid of the aligned plots
grid <- plot_grid(
  ncells_bar, logncells_bar, normncells_bar,
  nrow = 3,
  align = "v",
  axis = "tb"#,
  #labels = c("A", "B", "C")
)

# Plot ncells per cluster
grid

tiff(paste0("./figures/cell_type/var",r.variable,"sub_merged_cell_numbers_ratio_sorted_branch.tiff"),
     width = 30,height = 30,units = "cm", res = 300,compression = "lzw",bg = NA)
print(grid)

dev.off()




