# =============================================================================
# 05_cluster_selection_06_2025v3.R
# -----------------------------------------------------------------------------
# Purpose : Identify light/dark (LD) stress-responsive cell types. Summarises the
#           per-cell AUCell signature scores (IEGs, GRs, fkbp5) to the cell-type x
#           sample level, renders comparative heatmaps ordered by the cell-type
#           tree, and selects the top-quantile responsive clusters (Venn overlap)
#           to mask onto the UMAP. Ends by sourcing the neuropeptidergic subset.
# Inputs  : ./data/rds/step4_var<r.variable>.rds  (object with AUC scores)
#           ./data/rda/hcluster.rda (taxa_order); cell_type_table_m_modified.csv
# Outputs : AUC summary CSVs under ./outputs/; heatmap/Venn/UMAP PDFs under
#           ./figures/; updated `ss1` masks; sources 05_a for the nps subset.
# =============================================================================

## Load libraries

library(plyr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(harmony)
library(colorspace)
library(ggplot2)
library(viridis)
library(SCP)
#set the data and annotations
set.seed(1234)
options(future.globals.maxSize = 120000 * 1024^2)

r.variable=4000
vs="06_2025_v4"
umap="wnn.umap"

#set working dir (skip if already inside the version folder)
if(getwd()!=paste0("D:/projects/scMultiome_oxt/",vs)){
  setwd(paste0("./",vs))
}

#functions adopted and modified from Alberto Perez-Posada @apposada
source("./ext_code/r_code/functions/sourcefolder.R")
sourceFolder("./ext_code/r_code/functions/")

#other temp function: negated %in%
'%!in%' <- function(x,y)!('%in%'(x,y))
#load data
ss1=readRDS(paste0("./data/rds/step4_var",r.variable,".rds"))
load("./data/rda/hcluster.rda")
#AUC heatmap: pull the per-cell AUC scores + labels
libs=unique(ss1$orig.ident)
auc_matrix_tb= ss1[[c("merged_sub.anno_type","merged_sub_numb","orig.ident","AUC_ergs","AUC_GRs","AUC_fkbp5")]]

#auc : mean signature score per sample x cell type
auc_sum_matrix_tb=auc_matrix_tb %>%
  group_by(orig.ident, merged_sub_numb) %>%
  dplyr::summarise(mean_ergs=mean(AUC_ergs),
                   mean_GR=mean(AUC_GRs),
                   mean_fkbp5=mean(AUC_fkbp5))
#ucell_smooth (alternative: use UCell smoothed scores instead of AUCell)
'
auc_matrix_tb= ss1[[c("merged_sub.anno_type","orig.ident",
                      "ergsUcell_score_smooth_ucell",
                      "GRsUcell_score_smooth_ucell",
                      "fkbp5Ucell_score_smooth_ucell")]]
auc_sum_matrix_tb=auc_matrix_tb %>%
  group_by(orig.ident, merged_sub.anno_type) %>%
  dplyr::summarise(mean_ergs=mean(ergsUcell_score_smooth_ucell),
                   mean_GR=mean(GRsUcell_score_smooth_ucell),
                   mean_fkbp5=mean(fkbp5Ucell_score_smooth_ucell))
'
# Reshape to a cell-type x (signature x library) matrix, ordered by the tree
auc_sum_matrix=cbind(auc_sum_matrix_tb[c(which(auc_sum_matrix_tb$orig.ident==libs[1])),3:5],
                     auc_sum_matrix_tb[c(which(auc_sum_matrix_tb$orig.ident==libs[2])),3:5],
                     auc_sum_matrix_tb[c(which(auc_sum_matrix_tb$orig.ident==libs[3])),3:5],
                     auc_sum_matrix_tb[c(which(auc_sum_matrix_tb$orig.ident==libs[4])),3:5])
colnames(auc_sum_matrix)=paste0(c("AUC_ergs","AUC_GRs","AUC_fkbp5"),"_",rep(libs,each=3))
rownames(auc_sum_matrix)=unique(auc_sum_matrix_tb$merged_sub_numb)

texa_order_num=sapply(str_split(taxa_order[- grep("UnD",taxa_order)],"_"), function (x){x[[1]]})
texa_order_num=gsub("\\.0","",texa_order_num)
auc_sum_matrix=auc_sum_matrix[texa_order_num,]

#reorder into one matrix per signature (columns = 4 libraries)
auc_sum_matrix_erg=auc_sum_matrix[,c(1,4,7,10)]
auc_sum_matrix_GR=auc_sum_matrix[,c(2,5,8,11)]
auc_sum_matrix_fkbp5=auc_sum_matrix[,c(3,6,9,12)]

#visualization
#auc_heatmaps
##### IEG signature heatmap (cell types x libraries)

mt=as.matrix(auc_sum_matrix_erg)
clu_ha = HeatmapAnnotation(
  name = "library",
  cluster = factor(colnames(mt), levels = unique(colnames(mt))),
  col = list(cluster = setNames(rev(c("#d53e4f","#fee08b","#abdda4","#3288bd")) ,unique(colnames(mt)))),
  show_legend = F, show_annotation_name = F
)

auc_ergs <- ComplexHeatmap::Heatmap(
  name = "Auc_ergs",
  mt, #smed_wg_module_viz[,-modulecolumn],
  cluster_rows= F,
  show_row_names = T,
  row_names_side = "left",
  show_row_dend = F,
  cluster_columns = F,
  column_title = "IEGs",
  column_title_side = "bottom",
  show_column_names = F,
  column_names_side = "bottom",
  row_title_side = "left",
  row_title_rot = 0,
  border = T,border_gp = gpar(col = "grey", lty = 1),
  top_annotation = clu_ha,
  #bottom_annotation = clu_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = "grey", fill = NA, lty = 1))
  },
  use_raster = FALSE,
  col = colorRampPalette(c("white","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026"))(6)
)
draw(auc_ergs)


#auc_gr
##### GR signature heatmap

mt=as.matrix(auc_sum_matrix_GR)
clu_ha = HeatmapAnnotation(
  name = "library",
  cluster = factor(colnames(mt), levels = unique(colnames(mt))),
  col = list(cluster = setNames(rev(c("#d53e4f","#fee08b","#abdda4","#3288bd")) ,unique(colnames(mt)))),
  show_legend = F, show_annotation_name = F
)

auc_gr <- ComplexHeatmap::Heatmap(
  name = "Auc_GR",
  mt,
  cluster_rows= F,
  show_row_names = T,
  row_names_side = "left",
  show_row_dend = F,
  cluster_columns = F,
  column_title = "GRs",
  column_title_side = "bottom",
  show_column_names = F,
  column_names_side = "bottom",
  row_title_side = "left",
  row_title_rot = 0,
  border = T,border_gp = gpar(col = "grey", lty = 1),
  top_annotation = clu_ha,
  #bottom_annotation = clu_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = "grey", fill = NA, lty = 1))
  },
  use_raster = FALSE,
  col = colorRampPalette(c("white","#ffffcc","#d9f0a3","#addd8e","#238443","#005a32"))(6)
)
draw(auc_gr)


#auc_fkbp5
##### fkbp5 signature heatmap
mt=as.matrix(auc_sum_matrix_fkbp5)
clu_ha = HeatmapAnnotation(
  name = "library",
  cluster = factor(colnames(mt), levels = unique(colnames(mt))),
  col = list(cluster = setNames(rev(c("#d53e4f","#fee08b","#abdda4","#3288bd")) ,unique(colnames(mt)))),
  show_legend = F, show_annotation_name = F
)
auc_fkbp5 <- ComplexHeatmap::Heatmap(
  name = "Auc_fkbp5",
  mt,
  cluster_rows= F,
  show_row_names = T,
  row_names_side = "left",
  show_row_dend = F,
  cluster_columns = F,
  column_title = "fkbp5",
  column_title_side = "bottom",
  show_column_names = F,
  column_names_side = "bottom",
  row_title_side = "left",
  row_title_rot = 0,
  border = T,border_gp = gpar(col = "grey", lty = 1),
  top_annotation = clu_ha,
  #bottom_annotation = clu_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = "grey", fill = NA, lty = 1))
  },
  use_raster = FALSE,
  col = colorRampPalette(c("white","#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494"))(6)
)
draw(auc_fkbp5)



##### Combined three-signature heatmap panel

pdf(paste0("./figures/AUCh_heatmap.pdf"),
    width = 10,
    height = 30)
draw(auc_ergs+auc_gr+auc_fkbp5)

dev.off()

#all together (single clustered heatmap, split by signature)


mt=as.matrix(auc_sum_matrix)

clu_ha = HeatmapAnnotation(
  name = "library",
  cluster = factor(colnames(mt), levels = unique(colnames(mt))),
  col = list(cluster = setNames(c(rep("#d53e4f",3),rep("#fee08b",3),rep("#abdda4",3),rep("#3288bd",3)) ,unique(colnames(mt)))),
  show_legend = F, show_annotation_name = F
)


auc_all <- ComplexHeatmap::Heatmap(
  name = "Auc_all3",
  mt, #smed_wg_module_viz[,-modulecolumn],
  cluster_rows= T,clustering_method_rows = "ward.D2",
  show_row_names = T,row_names_side = "right",
  show_row_dend = T,row_dend_side = "left",
  cluster_columns = F,
  show_column_dend = F,
  show_column_names = TRUE,
  column_names_side = "bottom",
  row_title_side = "left",
  row_title_rot = 0,
  border = T,border_gp = gpar(col = "grey", lty = 1),
  column_split = rep(c("IEGs","GRs","fkbp5"),4),
  top_annotation = clu_ha,
  #bottom_annotation = clu_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = "grey", fill = NA, lty = 0.5))
  },
  use_raster = FALSE,row_km = 5,
  col = viridis(10)
)
a=draw(auc_all)

pdf(paste0("./figures/AUCh_all_heatmap_split.pdf"),
    width = 5,
    height = 12)
a

dev.off()



'
pdf(paste0("./figures/ucell_smooth_heatmap.pdf"),
     width = 10,
     height = 30)
draw(auc_ergs+auc_gr+auc_fkbp5)
dev.off()
'

### Summary statistics tables (mean +/- sd per cell type, with/without sample)

library(dplyr)
auc_stat_matrix_tb=as.tibble(auc_matrix_tb) %>%
  group_by(orig.ident, merged_sub.anno_type) %>%
  dplyr::summarise(n_total=n(),mean_ergs=mean(AUC_ergs),sd_ergs=sd(AUC_ergs),
                   mean_GR=mean(AUC_GRs),sd_GR=sd(AUC_GRs),
                   mean_fkbp5=mean(AUC_fkbp5),sd_fkbp5=sd(AUC_fkbp5))
write.csv(auc_stat_matrix_tb,"./outputs/auc_stat_matrix_tb.csv")

auc_stat_matrix_tb_sum=as.tibble(auc_matrix_tb) %>%
  group_by(merged_sub.anno_type) %>%
  dplyr::summarise(n_total=n(),mean_ergs=mean(AUC_ergs),sd_ergs=sd(AUC_ergs),
                   mean_GR=mean(AUC_GRs),sd_GR=sd(AUC_GRs),
                   mean_fkbp5=mean(AUC_fkbp5),sd_fkbp5=sd(AUC_fkbp5))
write.csv(auc_stat_matrix_tb_sum,"./outputs/auc_stat_matrix_tb_cluster_only.csv")

# Per-signature scaled means for the distribution scatter plot below
auc_stat_matrix_tb1=as.tibble(auc_matrix_tb) %>%
  group_by(merged_sub.anno_type,orig.ident) %>%
  dplyr::summarise(n_total=n(),mean_auc=mean(AUC_ergs))
auc_stat_matrix_tb1$mean_auc=scale(auc_stat_matrix_tb1$mean_auc)
auc_stat_matrix_tb1$type="IEGs"

auc_stat_matrix_tb2=as.tibble(auc_matrix_tb) %>%
  group_by(merged_sub.anno_type,orig.ident) %>%
  dplyr::summarise(n_total=n(),mean_auc=mean(AUC_GRs))
auc_stat_matrix_tb2$mean_auc=scale(auc_stat_matrix_tb2$mean_auc)
auc_stat_matrix_tb2$type="GRs"

auc_stat_matrix_tb3=as.tibble(auc_matrix_tb) %>%
  group_by(merged_sub.anno_type,orig.ident) %>%
  dplyr::summarise(n_total=n(),mean_auc=mean(AUC_fkbp5))
auc_stat_matrix_tb3$mean_auc=scale(auc_stat_matrix_tb3$mean_auc)
auc_stat_matrix_tb3$type="fkbp5"

auc_stat_matrix_tbs=rbind(auc_stat_matrix_tb1,auc_stat_matrix_tb2,auc_stat_matrix_tb3)


# Scatter: cluster size vs mean signal, per signature and sample
library(ggrepel)
dis_plot=ggplot(auc_stat_matrix_tbs, aes(x=log10(n_total), y=mean_auc,size=mean_auc,
                                         color= type,fill=type, label =merged_sub.anno_type,
                                         shape=orig.ident))
a=dis_plot+geom_point()+theme_classic()+ylim(c(0,6.5))+ geom_text_repel(size=3,nudge_y = 0.2)+ #geom_text(nudge_y = 0.25)+ #
  scale_shape_manual(values = c(21:24)) +
  scale_color_manual(values = c("fkbp5"="#2c7fb8","GRs"="#238443","IEGs"="#f03b20"))+
  scale_fill_manual(values = c("fkbp5"="#2c7fb8","GRs"="#238443","IEGs"="#f03b20"))


pdf(paste0("./figures/LD_responsive_distribution.pdf"),
    width = 8,height = 8)

print(a)
dev.off()

##### Select top-quantile responsive cell types per signature

#top quantile (cell types above the upper quartile for each signature)

top_q_list=list()
top_q_list["IEGs"]=auc_stat_matrix_tb_sum %>% dplyr::arrange(desc(mean_ergs))%>%
  dplyr::filter( mean_ergs> quantile(auc_stat_matrix_tb_sum$mean_ergs)[4])%>% dplyr::select(merged_sub.anno_type )
top_q_list["GRs"]=auc_stat_matrix_tb_sum %>% dplyr::arrange(desc(mean_GR))%>%
  dplyr::filter( mean_GR> quantile(auc_stat_matrix_tb_sum$mean_GR)[4])%>% dplyr::select(merged_sub.anno_type )
top_q_list["fkbp5"]=auc_stat_matrix_tb_sum %>% dplyr::arrange(desc(mean_fkbp5))%>%
  dplyr::filter(mean_fkbp5> quantile(auc_stat_matrix_tb_sum$mean_fkbp5)[4])%>% dplyr::select(merged_sub.anno_type )



###venndiagram : overlap of the three top-responsive sets
# Libraries
library(tidyverse)
library(hrbrthemes)
'%!in%' <- function(x,y)!('%in%'(x,y))

# library
library(ggvenn)

a=ggvenn(top_q_list, show_elements = T, stroke_size = 0.5,stroke_alpha = 0.5,
         label_sep = "\n", text_size = 3,fill_alpha = 0.1,
         fill_color =  c("#4401541A","#21908D1A","#FDE7251A" )
)
pdf(paste0("./figures/LD_responsive_venn.pdf"),
    width = 10,height = 10)

print(a)
dev.off()


#selection : keep the top responsive cell types per signature and mask on UMAP
IEG_top=top_q_list[[1]][c(1:6)]
GR_top=top_q_list[[2]][c(1:8)]
fkbp5_top=top_q_list[[3]][c(1:9)]

criteria.layers=c("IEG_top","GR_top","fkbp5_top")
mask_col=c("#f03b20","#238443","#2c7fb8")
for (i in 1:3) {
  top.l=get(criteria.layers[i])
  top_mask=as.character(ss1$merged_sub.anno_type)
  top_mask[c(which(ss1$merged_sub.anno_type %!in% top.l ))]=""   # blank out non-top cells
  ss1[[criteria.layers[i]]]=factor(top_mask,levels = unique(sort(top_mask))[-1])

  #umaps: highlight the selected cell types (with and without legend)
  ldr_map=CellDimPlot(
    srt = ss1, group.by = criteria.layers[i], seed = 0,
    reduction = umap, theme_use = "theme_blank",#legend.position = "none",
    palcolor = c(mask_col[i]),pt.size =0.7,
    label = T,label_insitu = T,label_repel = T
  )
  ldr_map2=CellDimPlot(
    srt = ss1, group.by = criteria.layers[i], seed = 0,
    reduction = umap, theme_use = "theme_blank",legend.position = "none",
    palcolor = c(mask_col[i]),pt.size = 0.7,
    label = F,label_insitu = T,label_repel = T
  )

  pdf(paste0("./figures/LD_responsive_mask_",criteria.layers[i],r.variable,".pdf"),
      width = 40,height = 20)

  print(ldr_map|ldr_map2)
  dev.off()
  pdf(paste0("./figures/LD_responsive_mask_2",criteria.layers[i],r.variable,".pdf"),
      width = 10,height = 10)

  print(ldr_map2)
  dev.off()
}

#####

#### Full annotated cell-type UMAPs (combined + per-sample)
type_table_m.v2= read.csv(paste0("./outputs/cell_type/cell_type_table_m_modified.csv"), header = T, row.names = 1)
type_table_m.v3= type_table_m.v2[- c(grep(pattern = "UnD", type_table_m.v2$numbered)),]
cycol=type_table_m.v2$clcol
names(cycol)=type_table_m.v2$numbered
a=DimPlot(ss1,group.by = "merged_sub.anno_type", seed = 0,split.by = "orig.ident",
          reduction = umap,cols = cycol,alpha = 1,pt.size = 0.6,
          label = F, label.size =3,stroke.size = 0.5,
          order=T
)+ ggtitle(NULL) & NoLegend()

b=DimPlot(ss1,group.by = "merged_sub.anno_type", seed = 0,
          reduction = umap,cols = cycol,alpha = 0.8,
          label = F, label.size =3.5,stroke.size = 0.5,#label.box = T,label.color = "white",
          order=T
)+ ggtitle(NULL) & NoLegend()

b=LabelClusters(b, id = "merged_sub.anno_type", repel = F,  fontface = "bold", color = "black"
                #,position = "nearest"
)

library(patchwork)
pdf(paste0("./figures/merged_sub_splited_umap_",r.variable,".pdf"), width = 50, height = 10)

print(((b|a)+plot_layout(widths = c(1, 4))))

dev.off()


ss1$merged_sub.anno_type=factor(ss1$merged_sub.anno_type, levels = str_sort(unique(ss1$merged_sub.anno_type),numeric = T))
b=DimPlot(ss1,group.by = "merged_sub.anno_type", seed = 0,
          reduction = umap,cols = cycol,alpha = 0.8,pt.size = 0.7,
          label = F, label.size =3.5,stroke.size = 0.5,#label.box = T,label.color = "white",
          order=T
)+ ggtitle(NULL)
b=LabelClusters(b, id = "merged_sub.anno_type", repel = F,  fontface = "bold", color = "black"
                #,position = "nearest"
)

pdf(paste0("./figures/merged_sub_umap_w_legened",r.variable,".pdf"), width = 15, height = 10)

print(b)

dev.off()


##GRs : glucocorticoid receptor expression on the UMAP

pdf(paste0("./figures/GRs_featureplot.pdf"),
    width = 10,height = 5)
p=FeatureDimPlot(
  srt = ss1, features = c("nr3c1","nr3c2"),
  #lower_cutoff = 0.8,
  assay = "RNA",
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor = c("lightgray","#ffffcc","#d9f0a3","#addd8e","#238443","#005a32"),
  seed = 0, compare_features = F, label =F, label_repel = T,label_insitu = TRUE,
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",#title = "GRs",
  theme_args=theme(strip.text = element_text(size = 20,face = c("bold.italic")) )
)
print(p)

dev.off()
#######
#for major NeuroEndo cells
####neuroendo cluster : subset and re-embed neuropeptidergic neurons (step 05_a)
source("./codes/05_a_nps_sub_umaps.R")

####

