
#Author: Min K Choi, m.choi@exter.ac.uk
#DEG analysis

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc()

## Load libraries

library(plyr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(colorspace)
library(SCP)
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

## Load Data
ss1=readRDS(paste0("./data/rds/Step6_var",r.variable,".rds"))

libs=unique(ss1$orig.ident)
type_table_m.v2= read.csv(paste0("./outputs/cell_type/cell_type_table_m_modified.csv"), header = T, row.names = 1)
type_table_m.v3= type_table_m.v2[- c(grep(pattern = "UnD", type_table_m.v2$numbered)),]

## Getting cluster information for pseudobulk
#scRNA-seq:
DefaultAssay(ss1)<-"RNA"

#pseudobulk for condition x clusters
ss1$merged_sub.anno_type_ori=paste0(ss1$merged_sub.anno_type,":",ss1$orig.ident)

modifID=gsub("_",replacement = ".",as.character(ss1$merged_sub.anno_type_ori))
modifID=gsub("/",replacement = ".",modifID)
modifID=gsub(":",replacement = ".",modifID)
modifID=gsub("\\+",replacement = "",modifID)

ss1[["modifID_cc"]]=modifID

Idents(ss1)=ss1$modifID_cc

# pseudo_split ------------------------------------------------------------

ss_idents <- ss1$modifID_cc

# replicates:
set.seed(1234)
ss_pseudo_reps <- sample(c("R1","R2"),length(ss_idents), replace = TRUE)
ss_pseudo_reps <- factor(ss_pseudo_reps, levels = c("R1","R2"))

# exp
exps <- rep("", length(ss_idents))

# the pseudobulk
psbulk_psrep_ss <- 
  pseudobulk_cond_rep(
    x = ss1@assays$RNA@counts,
    identities = ss_idents, 
    conditions = exps, 
    replicates = ss_pseudo_reps
  )

#tidyup
sampletable <- psbulk_psrep_ss$sampletable
sampletable$condition <- NULL
sampletable$id_combined <- sub("__","_",sampletable$id_combined)
sampletable$ctype <- factor(sampletable$ctype, levels = unique(ss_idents))
sampletable$replicate <- factor(sampletable$replicate, levels = unique(ss_pseudo_reps))
dac_mat <- psbulk_psrep_ss$matrix


####list_contrasts
post_cyte=sort(unique(ss_idents)[grep("post",unique(ss_idents))])
pre_cyte=sort(unique(ss_idents)[grep("pre",unique(ss_idents))])


list_contrasts <-
  as.list(as.data.frame(t(data.frame(
    "ctype",post_cyte,pre_cyte
  ))))
  

names(list_contrasts) <- sapply(list_contrasts, function(x){paste(x[c(2,3)], collapse = "_")})

list_contrasts


####list_contrasts2
bPAC_cyte=sort(unique(ss_idents)[grep("bPAC",unique(ss_idents))])
cont_cyte=sort(unique(ss_idents)[grep("cont",unique(ss_idents))])


list_contrasts2 <-
  as.list(as.data.frame(t(data.frame(
    "ctype",bPAC_cyte,cont_cyte
  ))))

names(list_contrasts2) <- sapply(list_contrasts2, function(x){paste(x[c(2,3)], collapse = "_")})

list_contrasts2


#loop
#####
dir.create("./outputs/DEGs")

postLD_preLD <- list()
for(i in names(list_contrasts)){
  postLD_preLD[[i]] <-
    deseq_sc(
      m = dac_mat,
      d = sampletable,
      contrast_info = list_contrasts[[i]],
      filter_by = "padj",
      p_threshold = 0.1,
      cell = i,
      plot_results = FALSE
    )
  write.csv(postLD_preLD[[i]]$res, paste0("./outputs/DEGs/DEGs_postLD_preLD_",gsub(":","_",i),".csv"))
}
y_max <- 25
library(gridExtra)
library(ggrastr)
library(ggrepel)
list_volcanos <- list()
for(celltype in names(postLD_preLD)){
  df <- postLD_preLD[[celltype]]$res
  df$padj[is.na(df$padj)]=1
  df$diff <- "none"
  df$diff[df$padj < 0.1 & df$log2FoldChange < 0] <- "down"
  df$diff[df$padj < 0.1 & df$log2FoldChange > 0] <- "up"
  df$diff <- factor(df$diff, levels = c("down","up","none"))
  
  df$label=rownames(df)
  df$label[which(df$padj>0.1)]=""
  df$label[which(df$pvalue > 0.05)]=""
  df$label[which(abs(df$log2FoldChange) < 1)]=""
  df$label[-c(unique(c(which(df$padj %in% sort(df$padj)[c(1:10)]),
                       which(df$log2FoldChange %in% rev(sort(df$log2FoldChange,))[c(1:10)])
  )))]=""
  
  list_volcanos[[celltype]] <- df %>%
    ggplot(aes(x = log2FoldChange, y = -log10(pvalue), color = diff)) +
    geom_point_rast(size = .75) + #geom_text_repel(size=3,max.overlaps = 20,force = 2,seed = 0)+
    scale_color_manual(values = c("none" = "#DAE7F2", "down" = "#174fbc", "up" ="firebrick" )) +
    ylim(c(0,y_max))+
    guides(color = "none") + 
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75), axis.title.x=element_blank(), axis.title.y=element_blank())+
    ggtitle(label = paste0(str_split(celltype,pattern = "_")[[1]][1],
                           "-pre"))
}

pdf(paste0("./figures/ss_diff_volcano_wo_label.pdf"), width = 40, height = 40)
grid.arrange(grobs = list_volcanos, ncol = 10)
dev.off()

#with label
list_volcanos_ft <- list()
for(celltype in names(postLD_preLD)){
  df <- postLD_preLD[[celltype]]$res
  df$diff <- "none"
  df$diff[df$pvalue < 0.05 & df$log2FoldChange < 0] <- "down"
  df$diff[df$pvalue < 0.05 & df$log2FoldChange > 0] <- "up"
  df$diff <- factor(df$diff, levels = c("down","up","none"))
  df$label=rownames(df)
  #df$label[is.na(df$padj)]=""
  df$label[which(df$padj>0.05)]=""
  df$label[which(df$pvalue > 0.05)]=""
  df$label[which(abs(df$log2FoldChange) < 1)]=""
  df$label[-c(unique(c(which(df$padj %in% sort(df$padj)[c(1:10)]),
                       which(df$log2FoldChange %in% rev(sort(df$log2FoldChange,))[c(1:10)])
  )))]=""
  
  #if (length(which(df$padj<0.1))>2) {
    list_volcanos_ft[[celltype]] <- df %>%
      ggplot(aes(x = log2FoldChange, y = -log10(pvalue), color = diff,label =label)) +
      geom_point_rast(size = .75) + geom_text_repel(size=4,max.overlaps = 10,min.segment.length = 2,seed = 0)+ #geom_text(size=5,nudge_y = 0.5) / geom_text_repel(size=4,max.overlaps = 10,min.segment.length = 2,seed = 0)+
      scale_color_manual(values = c("none" = "#DAE7F2", "down" = "#174fbc", "up" ="firebrick" )) +
      geom_hline(yintercept = -log10(0.05), col = "gray")+
      ylim(c(0,y_max))+
      guides(color = "none") + 
      theme_classic() +
      theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75), axis.title.x=element_blank(), axis.title.y=element_blank())+
      ggtitle(label = paste0(str_split(celltype,pattern = "_")[[1]][1],
                             "-pre"))
  #}
}

pdf(paste0("./figures/ss_diff_volcano_w_labels.pdf"), width = 30, height = 40)
grid.arrange(grobs = list_volcanos_ft, ncol = 10)
dev.off()




####
#DEG distribution
#pseudobulk DEG result
postLD_preLD_DEG_cl=c()
postLD_preLD_DEG_sample=c()
postLD_preLD_DEG_deg=c()
postLD_preLD_DEG_cell_count=c()

for (i in names(postLD_preLD)) {
  cl=sub("\\.([^.]*)$", "\\1",str_split(gsub("[^0-9.-]", "_", str_split(i,"\\_")[[1]][1]), "_", n=2)[[1]][1])
  sample=paste(last(str_split(str_split(i,"\\_")[[1]][1],"\\.")[[1]],n = 2)[1],
               last(str_split(str_split(i,"\\_")[[1]][1],"\\.")[[1]],n = 2)[2],sep="_")
  df_res.n=nrow(dplyr::filter(postLD_preLD[[i]]$res,padj < 0.1))
  sample_size=sum(postLD_preLD[[i]]$dds$ncells)
  postLD_preLD_DEG_cl=c(postLD_preLD_DEG_cl,cl)
  postLD_preLD_DEG_sample=c(postLD_preLD_DEG_sample,sample)
  postLD_preLD_DEG_deg=c(postLD_preLD_DEG_deg,df_res.n)
  postLD_preLD_DEG_cell_count=c(postLD_preLD_DEG_cell_count,sample_size)
}

postLD_preLD_tb= data.frame("cluster"=postLD_preLD_DEG_cl,
                            "comparison"=postLD_preLD_DEG_sample,
                            "n.DEG"=postLD_preLD_DEG_deg,
                            "cell.n.sum"=postLD_preLD_DEG_cell_count)

postLD_preLD_tb$logCount=log(postLD_preLD_tb$cell.n.sum,10) / 
  log(tapply(postLD_preLD_tb$cell.n.sum, postLD_preLD_tb$cluster, sum)[postLD_preLD_tb$cluster],10)
postLD_preLD_tb$comparison=gsub("_pre", "_LD_DEGs", postLD_preLD_tb$comparison)
postLD_preLD_tb$deg.ratio=postLD_preLD_tb$n.DEG*postLD_preLD_tb$logCount
write.csv(postLD_preLD_tb,paste0("./outputs/postLD_preLD_DEG_dist.tb.csv"))

#type_table_m.v2= read.csv(paste0("./outputs/cell_type/cell_type_table_m_modified.csv"), header = T, row.names = 1)
#type_table_m.v3= type_table_m.v2[- c(grep(pattern = "UnD", type_table_m.v2$numbered)),]


ptcol=type_table_m.v3$clcol
names(ptcol)=type_table_m.v3$cl_numb

library(ggrepel)
dis_plot=ggplot(postLD_preLD_tb, aes(x=logCount, y=n.DEG, size=n.DEG, 
                                     color= cluster, shape=comparison, label = paste(cluster,comparison, sep="_")
                                     ))
a=dis_plot+geom_point()+theme_classic()+guides(color = "none")+scale_color_manual(values = ptcol)  +geom_text_repel(size=5)

pdf(paste0("./figures/LD_nDEGs_post-pre_distribution.pdf"),
     width = 8,height = 8)
print(a)
dev.off()

load(file=paste0("./data/rda/hcluster.rda"))


DEG_matrix_tb= postLD_preLD_tb[,c("cluster","comparison","deg.ratio")]

DEG_matrix=cbind(DEG_matrix_tb[c(which(DEG_matrix_tb$comparison=="cont_post")),3],
                 DEG_matrix_tb[c(which(DEG_matrix_tb$comparison=="bPAC_post")),3])


colnames(DEG_matrix)=c("cont_DEGs","bPAC_DEGs")
rownames(DEG_matrix)=sapply(str_split(unique(DEG_matrix_tb$cluster),"_"), function(x){x[[1]]})

texa_order_num=sapply(str_split(taxa_order[- grep("UnD",taxa_order)],"_"), function (x){x[[1]]})
#texa_order_num=gsub("\\.0","",texa_order_num)
#texa_order_num[26]="34.0"
DEG_matrix=DEG_matrix[texa_order_num,]

#visualization
#auc_heatmaps
#####

#####
mt=as.matrix(DEG_matrix)
clu_ha = HeatmapAnnotation(
  name = "comparison",
  cluster = factor(colnames(mt), levels = unique(colnames(mt))),
  col = list(cluster = setNames(rev(c("#ec407a","#ffcdd2")) ,unique(colnames(mt)))),
  show_legend = F, show_annotation_name = F
)
DEG_prepostLD <- ComplexHeatmap::Heatmap(
  name = "post-pre.DEG",
  mt, 
  cluster_rows= F,
  show_row_names = T,
  row_names_side = "left",
  show_row_dend = F,
  cluster_columns = F,
  column_title = "post-pre",
  column_title_side = "bottom",
  column_title_rot = 90,
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
  col = colorRampPalette(c("white","#fce4ec","#f48fb1","#f06292","#e91e63","#c2185b"))(6)
)
draw(DEG_prepostLD)


######

#loop for bPAC vs cont
bPAC_cont <- list()
for(i in names(list_contrasts2)){
  bPAC_cont[[i]] <-
    deseq_sc(
      m = dac_mat,
      d = sampletable,
      contrast_info = list_contrasts2[[i]],
      filter_by = "padj",
      p_threshold = 0.1,
      cell = i,
      plot_results = FALSE
    )
  write.csv(bPAC_cont[[i]]$res, paste0("./outputs/DEGs/DEGs_bPAC_cont_",gsub(":","_",i),".csv"))
  
}



####volcano
#max_value <- max(unlist(sapply(hnf_DGE_all_broad, function(x){-log(x$res$pvalue)})))
#y_max <- as.numeric(as.character(cut(floor(max_value+10**magn_order(max_value)), breaks = c(0,10,30,50,75,100), labels = c(10,30,50,75,100))))

y_max <- 20
library(gridExtra)
library(ggrastr)
library(ggrepel)

#without label
list_volcanos <- list()
for(celltype in names(bPAC_cont)){
  df <- bPAC_cont[[celltype]]$res
  
  df$diff <- "none"
  df$diff[df$padj < 0.1 & df$log2FoldChange < 0] <- "down"
  df$diff[df$padj < 0.1 & df$log2FoldChange > 0] <- "up"
  df$diff <- factor(df$diff, levels = c("down","up","none"))
  
  df$label=rownames(df)
  df$label[which(df$padj>0.1)]=""
  df$label[which(df$pvalue > 0.05)]=""
  df$label[which(abs(df$log2FoldChange) < 1)]=""
  df$label[-c(unique(c(which(df$padj %in% sort(df$padj)[c(1:10)]),
                       which(df$log2FoldChange %in% rev(sort(df$log2FoldChange,))[c(1:10)])
  )))]=""
  
  list_volcanos[[celltype]] <- df %>%
    ggplot(aes(x = log2FoldChange, y = -log10(pvalue), color = diff)) +
    geom_point_rast(size = .75) + #geom_text_repel(size=3,max.overlaps = 20,force = 2,seed = 0)+
    scale_color_manual(values = c("none" = "#DAE7F2", "down" = "#174fbc", "up" ="firebrick" )) +
    ylim(c(0,y_max))+
    guides(color = "none") + 
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75), axis.title.x=element_blank(), axis.title.y=element_blank())+
    ggtitle(label = paste0(str_split(celltype,pattern = "_")[[1]][1],
                           "-cont"))
  
}

pdf(paste0("./figures/ss_diff_bPAC_cont_volcano.pdf"), width = 30, height = 40)
grid.arrange(grobs = list_volcanos, ncol = 10)
dev.off()

#with label
list_volcanos_ft <- list()
for(celltype in names(bPAC_cont)){
  df <- bPAC_cont[[celltype]]$res
  df$diff <- "none"
  df$diff[df$padj < 0.1 & df$log2FoldChange < 0] <- "down"
  df$diff[df$padj < 0.1 & df$log2FoldChange > 0] <- "up"
  df$diff <- factor(df$diff, levels = c("down","up","none"))
  df$label=rownames(df)
  df$label[is.na(df$padj)]=""
  df$label[which(df$padj>0.05)]=""
  df$label[which(df$pvalue > 0.05)]=""
  df$label[which(abs(df$log2FoldChange) < 1)]=""
  df$label[-c(unique(c(which(df$padj %in% sort(df$padj)[c(1:10)]),
                       which(df$log2FoldChange %in% rev(sort(df$log2FoldChange,))[c(1:10)])
  )))]=""
  
  #if (length(which(df$padj<0.1))>2) {
    list_volcanos_ft[[celltype]] <- df %>%
      ggplot(aes(x = log2FoldChange, y = -log10(pvalue), color = diff,label =label)) +
      geom_point_rast(size = .75) + geom_text_repel(size=4,max.overlaps = 10,min.segment.length = 2,seed = 0)+ #geom_text(size=5,nudge_y = 0.5) / geom_text_repel(size=4,max.overlaps = 10,min.segment.length = 2,seed = 0)+
      scale_color_manual(values = c("none" = "#DAE7F2", "down" = "#174fbc", "up" ="firebrick" )) +
      geom_hline(yintercept = -log10(0.05), col = "gray")+
      ylim(c(0,y_max))+
      guides(color = "none") + 
      theme_classic() +
      theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75), axis.title.x=element_blank(), axis.title.y=element_blank())+
      ggtitle(label = paste0(str_split(celltype,pattern = "_")[[1]][1],
                             "-cont"))
  #}
}

pdf(paste0("./figures/ss_diff_bPAC_cont_volcano_with_label.pdf"), width = 30, height = 40)
grid.arrange(grobs = list_volcanos_ft,ncol=10)
dev.off()

####
#DEG distribution
#pseudobulk DEG result
bPAC_cont_DEG_cl=c()
bPAC_cont_DEG_sample=c()
bPAC_cont_DEG_deg=c()
bPAC_cont_DEG_cell_count=c()

for (i in names(bPAC_cont)) {
  cl=sub("\\.([^.]*)$", "\\1",str_split(gsub("[^0-9.-]", "_", str_split(i,"\\_")[[1]][1]), "_", n=2)[[1]][1])
  sample=paste(last(str_split(str_split(i,"\\_")[[1]][1],"\\.")[[1]],n = 2)[1],
               last(str_split(str_split(i,"\\_")[[1]][1],"\\.")[[1]],n = 2)[2],sep="_")
  df_res.n=nrow(dplyr::filter(bPAC_cont[[i]]$res,padj < 0.1))
  sample_size=sum(bPAC_cont[[i]]$dds$ncells)
  bPAC_cont_DEG_cl=c(bPAC_cont_DEG_cl,cl)
  bPAC_cont_DEG_sample=c(bPAC_cont_DEG_sample,sample)
  bPAC_cont_DEG_deg=c(bPAC_cont_DEG_deg,df_res.n)
  bPAC_cont_DEG_cell_count=c(bPAC_cont_DEG_cell_count,sample_size)
}


bPAC_cont_tb= data.frame("cluster"=bPAC_cont_DEG_cl,
                            "comparison"=bPAC_cont_DEG_sample,
                            "n.DEG"=bPAC_cont_DEG_deg,
                            "cell.n.sum"=bPAC_cont_DEG_cell_count)

bPAC_cont_tb$logCount=log(bPAC_cont_tb$cell.n.sum,10) / log(tapply(bPAC_cont_tb$cell.n.sum, bPAC_cont_tb$comparison, sum)[bPAC_cont_tb$comparison],10)
bPAC_cont_tb$deg.ratio=bPAC_cont_tb$n.DEG*bPAC_cont_tb$logCount

#bPAC_cont_tb$comparison=gsub("_pre", "_LD_DEGs", bPAC_cont_tb$comparison)
write.csv(bPAC_cont_tb,paste0("./outputs/bPAC_cont_DEG_dist.tb.csv"))



ptcol=type_table_m.v3$clcol
names(ptcol)=type_table_m.v3$cl_numb

library(ggrepel)
dis_plot=ggplot(bPAC_cont_tb, aes(x=logCount, y=n.DEG, size=n.DEG, 
                                     color= cluster, shape=comparison, label = paste(cluster,comparison, sep="_")
))
a=dis_plot+geom_point()+theme_classic()+guides(color = "none")+scale_color_manual(values = ptcol)  +geom_text_repel(size=5)

pdf(paste0("./figures/geno_nDEGs_distribution.pdf"),
    width = 8,height = 8)
print(a)
dev.off()

load(file=paste0("./data/rda/hcluster.rda"))


DEG_matrix_tb2= bPAC_cont_tb[,c("cluster","comparison","deg.ratio")]

DEG_matrix2=cbind(DEG_matrix_tb2[c(which(DEG_matrix_tb2$comparison=="bPAC_pre")),3],
                 DEG_matrix_tb2[c(which(DEG_matrix_tb2$comparison=="bPAC_post")),3])


colnames(DEG_matrix2)=c("bPAC_pre","bPAC_DEGs")
rownames(DEG_matrix2)=sapply(str_split(unique(DEG_matrix_tb2$cluster),"_"), function(x){x[[1]]})

texa_order_num=sapply(str_split(taxa_order[- grep("UnD",taxa_order)],"_"), function (x){x[[1]]})
#texa_order_num=gsub("\\.0","",texa_order_num)
#texa_order_num[26]="34.0"
DEG_matrix2=DEG_matrix2[texa_order_num,]

#visualization
#auc_heatmaps
#####

#####
mt=as.matrix(DEG_matrix2)
clu_ha = HeatmapAnnotation(
  name = "comparison",
  cluster = factor(colnames(mt), levels = unique(colnames(mt))),
  col = list(cluster = setNames(rev(c("#6a1b9a","#d1c4e9")) ,unique(colnames(mt)))),
  show_legend = F, show_annotation_name = F
)
DEG_bPACCont <- ComplexHeatmap::Heatmap(
  name = "bPAC-cont.DEG",
  mt, 
  cluster_rows= F,
  show_row_names = T,
  row_names_side = "left",
  show_row_dend = F,
  cluster_columns = F,
  column_title = "bPAC-contr.",
  column_title_side = "bottom",
  column_title_rot = 90,
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
  col = colorRampPalette(c("white","#e1bee7","#ba68c8","#9c27b0","#7b1fa2","#6a1b9a"))(6)
)
draw(DEG_bPACCont)

#####
auc_matrix_tb= ss1[[c("merged_sub_numb","orig.ident","AUC_ergs","AUC_GRs","AUC_fkbp5")]]

#auc
auc_sum_matrix_tb=auc_matrix_tb %>% 
  group_by(orig.ident, merged_sub_numb) %>% 
  dplyr::summarise(mean_ergs=mean(AUC_ergs),
                   mean_GR=mean(AUC_GRs),
                   mean_fkbp5=mean(AUC_fkbp5))
#ucell_smooth
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
auc_sum_matrix=cbind(auc_sum_matrix_tb[c(which(auc_sum_matrix_tb$orig.ident==libs[1])),3:5],
                     auc_sum_matrix_tb[c(which(auc_sum_matrix_tb$orig.ident==libs[2])),3:5],
                     auc_sum_matrix_tb[c(which(auc_sum_matrix_tb$orig.ident==libs[3])),3:5],
                     auc_sum_matrix_tb[c(which(auc_sum_matrix_tb$orig.ident==libs[4])),3:5])
colnames(auc_sum_matrix)=paste0(c("AUC_ergs","AUC_GRs","AUC_fkbp5"),"_",rep(libs,each=3))
rownames(auc_sum_matrix)=unique(auc_sum_matrix_tb$merged_sub_numb)

texa_order_num=sapply(str_split(taxa_order[- grep("UnD",taxa_order)],"_"), function (x){x[[1]]})
texa_order_num=gsub("\\.0","",texa_order_num)
auc_sum_matrix=auc_sum_matrix[texa_order_num,]

#reorder
auc_sum_matrix_erg=auc_sum_matrix[,c(1,4,7,10)]
auc_sum_matrix_GR=auc_sum_matrix[,c(2,5,8,11)]
auc_sum_matrix_fkbp5=auc_sum_matrix[,c(3,6,9,12)]

#visualization
#auc_heatmaps
#####

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
#####

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
#####
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
#####

pdf(paste0("./figures/AUC_DEG_heatmap.pdf"),
     width = 5,
     height = 15)
draw(auc_ergs+auc_gr+auc_fkbp5+DEG_prepostLD+DEG_bPACCont)
dev.off()

pdf(paste0("./figures/DEG_heatmap.pdf"),
     width = 4,
     height = 16)
draw(DEG_prepostLD+DEG_bPACCont)

dev.off()


save(
  list_contrasts,
  list_contrasts2,
  dac_mat,
  sampletable,
  postLD_preLD,
  bPAC_cont,
  file=paste0("./data/rda/DEGs_step7.rda")
)

#s.tables for DEG

postLD_preLD_DEG_All=postLD_preLD[[names(list_contrasts)[1]]]$res
postLD_preLD_DEG_All[,"comparison"]=names(list_contrasts)[1]
for (i in names(list_contrasts)[-1]) {
  tb=postLD_preLD[[i]]$res
  tb[,"comparison"]=i
  postLD_preLD_DEG_All=rbind(postLD_preLD_DEG_All,tb)
}

postLD_preLD_DEG_All_ft= postLD_preLD_DEG_All %>% dplyr::filter(padj <0.1) %>% dplyr::filter(abs(log2FoldChange) > log2(1.5))

bPAC_cont_DEG_All=bPAC_cont[[names(list_contrasts2)[1]]]$res
bPAC_cont_DEG_All[,"comparison"]=names(list_contrasts2)[1]

for (i in names(list_contrasts2)[-1]) {
  tb=bPAC_cont[[i]]$res
  tb[,"comparison"]=i
  bPAC_cont_DEG_All=rbind(bPAC_cont_DEG_All,tb)
}

bPAC_cont_DEG_All_ft= bPAC_cont_DEG_All %>% dplyr::filter(padj <0.1) %>% dplyr::filter(abs(log2FoldChange) > log2(1.5))

DEG_All_ft=rbind(postLD_preLD_DEG_All_ft,bPAC_cont_DEG_All_ft)
write.csv(DEG_All_ft,"./outputs/DEGs/DEG_ft_all.csv")
#####
#ATAC-DAR

DefaultAssay(ss1)="ATAC_macs3"
# the pseudobulk
psbulk_psrep_scatac <- 
  pseudobulk_cond_rep(
    x = ss1@assays$ATAC_macs3@counts,
    identities = ss_idents, 
    conditions = exps, 
    replicates = ss_pseudo_reps
  )


#tidyup
sampletable <- psbulk_psrep_scatac$sampletable
sampletable$condition <- NULL
sampletable$id_combined <- sub("__","_",sampletable$id_combined)
sampletable$ctype <- factor(sampletable$ctype, levels = unique(ss_idents))
sampletable$replicate <- factor(sampletable$replicate, levels = unique(ss_pseudo_reps))
dac_mat <- psbulk_psrep_scatac$matrix

#
####

list_contrasts

list_contrasts2


#loop
#####
postLD_preLD_atac <- list()
for(i in names(list_contrasts)){
  postLD_preLD_atac[[i]] <-
    deseq_sc(
      m = dac_mat,
      d = sampletable,
      contrast_info = list_contrasts[[i]],
      filter_by = "pvalue",
      p_threshold = 0.1,
      cell = i,
      plot_results = FALSE
    )
}

bPAC_cont_ATAC <- list()
for(i in names(list_contrasts2)){
  bPAC_cont_ATAC[[i]] <-
    deseq_sc(
      m = dac_mat,
      d = sampletable,
      contrast_info = list_contrasts2[[i]],
      filter_by = "pvalue",
      p_threshold = 0.1,
      cell = i,
      plot_results = FALSE
    )
}

gene_cordi=ss1@assays$ATAC_macs3@annotation
for(i in names(postLD_preLD_atac)){
  peak_coords <- postLD_preLD_atac[[i]]$diffgenes

  if(length(peak_coords)>1){
    ranges.show <- StringToGRanges(peak_coords)
    hits <- findOverlaps(gene_cordi,ranges.show)
    selected_elements <- gene_cordi[queryHits(hits)]
   postLD_preLD_atac[[i]]$diff_gene_link <-
    unique(selected_elements$gene_name)
  }else if(length(peak_coords)==1){
    if(is.na(peak_coords) == F){
      ranges.show <- StringToGRanges(peak_coords)
      hits <- findOverlaps(gene_cordi,ranges.show)
      selected_elements <- gene_cordi[queryHits(hits)]
      postLD_preLD_atac[[i]]$diff_gene_link <-
        unique(selected_elements$gene_name)
    }
  }
  print(i)
}


for(i in names(bPAC_cont_ATAC)){
  peak_coords <- bPAC_cont_ATAC[[i]]$diffgenes
  
  if(length(peak_coords)>1){
    ranges.show <- StringToGRanges(peak_coords)
    hits <- findOverlaps(gene_cordi,ranges.show)
    selected_elements <- gene_cordi[queryHits(hits)]
    bPAC_cont_ATAC[[i]]$diff_gene_link <-
      unique(selected_elements$gene_name)
  }else if(length(peak_coords)==1){
    if(is.na(peak_coords) == F){
      ranges.show <- StringToGRanges(peak_coords)
      hits <- findOverlaps(gene_cordi,ranges.show)
      selected_elements <- gene_cordi[queryHits(hits)]
      bPAC_cont_ATAC[[i]]$diff_gene_link <-
        unique(selected_elements$gene_name)
    }
  }
  print(i)
}

#DAR
library(gridExtra)
library(ggrastr)
library(ggrepel)

y_max <- 10
list_volcanos <- list()
for(celltype in names(postLD_preLD_atac)){
  df <- postLD_preLD_atac[[celltype]]$res
  
  df$diff <- "none"
  df$diff[df$pvalue < 0.05 & df$log2FoldChange < 0] <- "down"
  df$diff[df$pvalue < 0.05 & df$log2FoldChange > 0] <- "up"
  df$diff <- factor(df$diff, levels = c("down","up","none"))
  
  df$label=rownames(df)
  df$label[which(df$pvalue > 0.05)]=""
  df$label[which(abs(df$log2FoldChange) < 1)]=""
  df$label[-c(unique(c(which(df$pvalue %in% sort(df$pvalue)[c(1:10)]),
                       which(df$log2FoldChange %in% rev(sort(df$log2FoldChange,))[c(1:10)])
  )))]=""
  
  list_volcanos[[celltype]] <- df %>%
    ggplot(aes(x = log2FoldChange, y = -log10(pvalue), color = diff)) +
    geom_point_rast(size = .75) + #geom_text_repel(size=3,max.overlaps = 20,force = 2,seed = 0)+
    scale_color_manual(values = c("none" = "#DAE7F2", "down" = "#174fbc", "up" ="firebrick" )) +
    ylim(c(0,y_max))+
    guides(color = "none") + 
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75), axis.title.x=element_blank(), axis.title.y=element_blank())+
    ggtitle(label = paste0(str_split(celltype,pattern = "_")[[1]][1],
                           "-pre"))
}

pdf(paste0("./figures/ss_diff_LD_volcano_ATAC_wo_label.pdf"), width = 30, height = 40)
grid.arrange(grobs = list_volcanos, ncol = 10)
dev.off()

y_max=15
list_volcanos_geno <- list()
for(celltype in names(bPAC_cont_ATAC)){
  df <- bPAC_cont_ATAC[[celltype]]$res
  df$diff <- "none"
  df$diff[df$pvalue < 0.05 & df$log2FoldChange < 0] <- "down"
  df$diff[df$pvalue < 0.05 & df$log2FoldChange > 0] <- "up"
  df$diff <- factor(df$diff, levels = c("down","up","none"))
  df$label=rownames(df)
  df$label[is.na(df$pvalue)]=""
  df$label[which(df$pvalue > 0.05)]=""
  df$label[which(abs(df$log2FoldChange) < 1)]=""
  df$label[-c(unique(c(which(df$pvalue %in% sort(df$pvalue)[c(1:10)]),
                       which(df$log2FoldChange %in% rev(sort(df$log2FoldChange,))[c(1:10)])
  )))]=""
  
  if (length(which(df$pvalue< 0.05))>2) {
    list_volcanos_geno[[celltype]] <- df %>%
      ggplot(aes(x = log2FoldChange, y = -log10(pvalue), color = diff,label =label)) +
      geom_point_rast(size = .75) + #geom_text_repel(size=4,max.overlaps = 10,min.segment.length = 2,seed = 0)+ #geom_text(size=5,nudge_y = 0.5) / geom_text_repel(size=4,max.overlaps = 10,min.segment.length = 2,seed = 0)+
      scale_color_manual(values = c("none" = "#DAE7F2", "down" = "#174fbc", "up" ="firebrick" )) +
      geom_hline(yintercept = -log10(0.05), col = "gray")+
      ylim(c(0,y_max))+
      guides(color = "none") + 
      theme_classic() +
      theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.75), axis.title.x=element_blank(), axis.title.y=element_blank())+
      ggtitle(label = paste0(str_split(celltype,pattern = "_")[[1]][1],
                             "-cont"))
  }
}
pdf(paste0("./figures/ss_diff_geno_volcano_ATAC_wo_label.pdf"), width = 30, height = 40)
grid.arrange(grobs = list_volcanos_geno, ncol = 10)
dev.off()


save(
  list_contrasts,
  list_contrasts2,
  dac_mat,
  sampletable,
  postLD_preLD_atac,
  bPAC_cont_ATAC,
file=paste0("./data/rda/DAR_ATAC.rda")
  )

saveRDS(ss1,paste0("./data/rds/Step7_var",r.variable,".rds"))

#s.tables for DEG

postLD_preLD_DAR_All=postLD_preLD_atac[[names(list_contrasts)[1]]]$res
postLD_preLD_DAR_All[,"comparison"]=names(list_contrasts)[1]
for (i in names(list_contrasts)[-1]) {
  tb=postLD_preLD_atac[[i]]$res
  tb[,"comparison"]=i
  postLD_preLD_DAR_All=rbind(postLD_preLD_DAR_All,tb)
}

postLD_preLD_DAR_All_ft= postLD_preLD_DAR_All %>% dplyr::filter(pvalue <0.05) %>% dplyr::filter(abs(log2FoldChange) > log2(1.2))

bPAC_cont_DAR_All=bPAC_cont_ATAC[[names(list_contrasts2)[1]]]$res
bPAC_cont_DAR_All[,"comparison"]=names(list_contrasts2)[1]
for (i in names(list_contrasts2)[-1]) {
  tb=bPAC_cont_ATAC[[i]]$res
  tb[,"comparison"]=i
  bPAC_cont_DAR_All=rbind(bPAC_cont_DAR_All,tb)
}

bPAC_cont_DAR_All_ft= bPAC_cont_DAR_All %>% dplyr::filter(pvalue <0.05) %>% dplyr::filter(abs(log2FoldChange) > log2(1.2))

DAR_All_ft=rbind(postLD_preLD_DAR_All_ft,bPAC_cont_DAR_All_ft)
write.csv(DAR_All_ft,"./outputs/DEGs/DAR_ft_all.csv")


######

#####DEG_OCR for cell_group=c("35.0_avp.crhb", "35.1_avp","45_sst1.1")

#load(file=paste0("./data/rda/DEGs_step7.rda"))

library(presto)
library(GenomicRanges)
library(viridis)
library(Signac)

DefaultAssay(ss1) <- "ATAC_macs3"
dir.create("./figures/ATAC/DEGs/")

#cell_group=c("35.0_avp.crhb", "35.1_avp","45_sst1.1")
cell_group=c("35.0", "35.1","45")
for (i in cell_group) {
  
  target_groups=unique(
    c(names(postLD_preLD)[grepl(pattern = i,names(postLD_preLD))], names(bPAC_cont)[grepl(pattern = i,names(bPAC_cont))])
    )
  
  l.deg=c()
  
  for(c.deg in target_groups) {
    if (c.deg %in% names(postLD_preLD)) {
      l.deg=c(l.deg,unique(c(postLD_preLD[[c.deg]]$diffgene)))
    }else{
      l.deg=c(l.deg,unique(c(bPAC_cont[[c.deg]]$diffgene)))
    }
  }
  l.deg=na.omit(l.deg)
  subgroup2show=unique(ss1$merged_sub.anno_type[grepl(pattern = i,ss1$merged_sub.anno_type)])
  
  #top dar
  DA_ct_npc= wilcoxauc(ss1,
                       groups_use = unique(ss1$modifID_cc[grepl(pattern = i,ss1$modifID_cc)]),
                       group_by = "modifID_cc", 
                       seurat_assay = "ATAC_macs3")
  
  top_peaks_ct <- DA_ct_npc %>%
    dplyr::filter(abs(logFC) > log2(1.25) &
                    padj < 0.1 
                  #pct_in - pct_out > 13 &
                  #auc > 0.55
    ) %>%
    group_by(group) 
  
  ranges.show <- StringToGRanges(top_peaks_ct$feature)
  ranges.show$color <- "gray"
  
  
  for (k in l.deg) {
    gene_cordi=LookupGeneCoords(ss1,k)
    hits <- findOverlaps(ranges.show, gene_cordi)
    
    # Extract elements that fall within the range
    selected_elements <- ranges.show[queryHits(hits)]
    
    a=CoveragePlot(
      object = ss1,
      split.by = "orig.ident",#peaks.group.by = "orig.ident",  
      group.by = "merged_sub.anno_type",  
      region = gene_cordi,
      features = k,
      region.highlight = selected_elements,
      links = T,
      #heights = c(18,1,1,1),
      expression.assay = "RNA",
      idents = subgroup2show,
      extend.upstream = 2000,
      extend.downstream = 2000
    )
    
    
    n=str_split(str_split(string = a[[1]][[1]][["labels"]][["y"]],pattern = "- ")[[1]][2],"\\)")[[1]][1]
    if (as.integer(n) > 4) {
      pdf(paste0("./figures/ATAC/DEGs/",i,"_ATAC_Degs_",k,".pdf"),
           width = 7.5,height = 10)
      print(a & scale_fill_manual(values = magma(5,alpha = 0.5)))
      dev.off()
      aa=paste0(k," has enough peak")
      print(aa)
    }else{
      aa=paste0(k," not enough peak")
      print(aa)
    }
  }
}

