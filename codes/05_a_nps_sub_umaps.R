
ss2=subset(ss1, np100_clusters_mask_simple != "_" )
DefaultAssay(ss2)="RNA"
DimPlot(ss2,group.by = "merged_sub_numb", seed = 0,
        reduction = umap,alpha = 0.8,pt.size = 0.7,
        label = T, label.size =3,stroke.size = 0.5,
        order=T
)+ ggtitle(NULL)& NoLegend()


ss2 <- FindMultiModalNeighbors(ss2, reduction.list = list("pca", "lsi"), 
                                          dims.list = list(1:100, 2:30))
ss2 <- RunUMAP(ss2, nn.name = "weighted.nn", n.neighbors = 45,seed.use =362,
                          metric = "cosine",n.components = 2, min.dist = 0.2,alpha = 1, gamma = 1.0,
                          reduction.name = "wnn.umap2", reduction.key = "wnnUMAP2_")

library(ggrepel)

ctcol=unique(ss2$ctype_col)
names(ctcol)=unique(ss2$merged_sub.anno_type)
ss2$merged_sub.anno_type=factor(ss2$merged_sub.anno_type, levels = str_sort(unique(ss2$merged_sub.anno_type), numeric = T))
p7=DimPlot(ss2, reduction = "wnn.umap2", group.by = "merged_sub.anno_type", 
        cols=ctcol,#label.box = T,label.color = "white",
        label = F, repel = TRUE)+ ggtitle(NULL)& NoLegend()

p7.1=LabelClusters(p7, id = "merged_sub.anno_type", repel = T,  
              fontface = "bold", color = "darkgrey",
              nudge_x = 2,nudge_y = -1,
              box = F#,
              #label.padding=0.5,label.r=0.5
              #position = "nearest"
)
p7.1
pdf(paste0("./figures/cell_type/var",r.variable,"_ss2_np50_mask_umap.pdf"),
     width = 5,height = 5)
print(p7.1)

dev.off()

nps=c("avp","crhb","sst1.1","oxt","galn","prlh2","trh","nts","th","th2","kiss1","npy")

p8=FeaturePlot(ss2,features = nps,ncol = 3,
            reduction = "wnn.umap2",pt.size = 0.7,order = T,min.cutoff = 0.8
)

pdf(paste0("./figures/cell_type/var",r.variable,"_ss2_np50_feature_umap.pdf"),
     width = 12,height = 15)
print(p8)

dev.off()

ctcol=ctcol[str_sort(names(ctcol),numeric = T)]
scp_p0=FeatureStatPlot(ss2, stat.by = nps,
                       assay = "RNA",#split.by = "orig.ident",
                       add_box = T,
                       #bg.by = "merged_sub.anno",bg_alpha = 0.1,
                       #add_stat = c("median"),#c("none", "mean", "median"),
                       stat_color = "grey",stat_size = 0.5,
                       ylab = "",xlab = "",#bg_palette =  "simspec",
                       palcolor = ctcol,
                       sort = T,legend.position = "none",
                       add_line = 0.5,line_size = 0.5,#line_color = "gray",
                       #stat_stroke = 0.5,stat_shape = 3,
                       #legend.position = "top", legend.direction = "horizontal",
                       group.by = "merged_sub_numb", 
                       plot.by ="group",   #"feature", 
                       stack =T)

pdf(paste("./figures/cell_type/nps_mk_r4000_scp.pdf"),
     width = 10,height = 10)
print(scp_p0)
dev.off()

###
ss2=PrepSCTFindMarkers(ss2)
#identifiy markers
ss2.markers= FindAllMarkers(
  ss2, assay = "SCT",
  group.by = "merged_sub.anno_type",
  logfc.threshold = 0.1,
  test.use = "MAST",#"DESeq2","MAST","wilcox"
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
top100_mk= ss2.markers %>% 
  dplyr::filter(abs(avg_log2FC)>1) %>%
  dplyr::filter(p_val_adj <0.05) %>%
  group_by(cluster) %>%
  top_n(100, wt=abs(avg_log2FC))

write.csv(ss2.markers,paste0("./outputs/cell_type/nps.all.merged.markers_",
                             "var",r.variable,
                             ".csv"))

saveRDS(ss2, "./data/rds/step5_nps_var4000.rds")
####cont_pre only
ss2.ori=ss2

##markers among nps
ss2=subset(ss2.ori, orig.ident == "cont_pre")
ss2=PrepSCTFindMarkers(ss2)
#identifiy markers
ss2.markers= FindAllMarkers(
  ss2, assay = "SCT",
  group.by = "merged_sub.anno_type",
  logfc.threshold = 0.1,
  test.use = "MAST",#"DESeq2","MAST","wilcox"
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
top100_mk= ss2.markers %>% 
  dplyr::filter(abs(avg_log2FC)>1) %>%
  dplyr::filter(p_val_adj <0.05) %>%
  group_by(cluster) %>%
  top_n(100, wt=abs(avg_log2FC))

write.csv(ss2.markers,paste0("./outputs/cell_type/ss2.merged.markers_",
                                    "var",r.variable,
                                    ".csv"))
#ss2.markers=read.csv(paste0("./outputs/cell_type/ss2.merged.markers_","var",r.variable,".csv"))

#top50s
top_mk.hc= ss2.markers %>% 
  dplyr::filter(abs(avg_log2FC)>1) %>%
  dplyr::filter(p_val_adj <0.1) %>%
  dplyr::filter(pct.1 >0.3) %>%
  group_by(cluster) %>%
  top_n(20, wt= -p_val_adj )

top_mk.hc$cluster=factor(top_mk.hc$cluster, levels = str_sort(unique(top_mk.hc$cluster), numeric = T))


mk2_list=list()
dir.create(paste0("./figures/cell_type/nps_mks"))
for (i in 1:length(unique(top_mk.hc$cluster))) {
  p=FeatureDimPlot(
    srt = ss2.ori, features = top_mk.hc$gene[which(top_mk.hc$cluster == unique(top_mk.hc$cluster)[i])],
    assay = "SCT",
    seed = 0, compare_features = F, label =F, label_repel = T,label_insitu = TRUE, 
    add_density = F,palette = "GdRd",bg_cutoff = 0.5, 
    reduction = "wnn.umap2", 
    theme_use = "theme_blank",
    theme_args=theme(strip.text = element_text(size = 20,face = c("bold.italic")))
  )
  mk2_list[[i]]=p+ #ggtitle(label ="", subtitle = paste0("cl.=",cluster.nm[i]))+
    theme(title = element_text(size = 0),
          plot.subtitle=element_text(size = 10))
  
  pdf(paste0("./figures/cell_type/nps_mks/nps_cont_pre_ss_feature_umap_",unique(top_mk.hc$cluster)[i],".pdf"),
      width = 25,height = 20)
  plot(p)
  
  dev.off()
}


npcols=unique(ss2$ctype_col)
names(npcols)=unique(ss2$merged_sub.anno_type)
npcols=npcols[str_sort(unique(top_mk.hc$cluster), numeric = T)]

npc=as.character(npcols)

#top10
levels(ss2$merged_sub.anno_type)=

ht5 <- GroupHeatmap(ss2,assay = "RNA",
                    features = top_mk.hc$gene, feature_split = top_mk.hc$cluster, 
                    group.by = "merged_sub.anno_type",#split.by = "orig.ident", 
                    heatmap_palette = "YlOrRd",
                    #group_palcolor ="",# "#FFD9A4" "#FF8F6C" "#EEBDD1" "#9DE7D7" "#FFE671"
                    feature_split_palcolor=as.list(npcols),
                    #cell_split_palcolor = rev(c("#d53e4f","#fee08b","#abdda4","#3288bd")),
                    group_palcolor = list(c(npc)),column_title = "NP Cell clusters",                    
                    #group_palette = "viridis",
                    #cell_annotation = c(nps),#,"Phase","foxp2","foxp4", "oxt"), 
                    #cell_annotation_palette = simcol,
                    #cell_annotation_params = list(height = unit(10, "mm")),
                    feature_annotation = c("TF","dGCprimed"),label_size = 1.5,#flip = T,
                    feature_annotation_palcolor = list(c("gold", "steelblue"),c("red")),
                    add_dot = T, add_bg = TRUE, nlabel = 0, show_row_names = TRUE
)

print(ht5$plot)

pdf(paste0("./figures/cell_type/hcluster_top10nps_orig_mk_r4000_mergedsub_wlabel.pdf"),
     width = 15,height = 50)
print(ht5$plot)
dev.off()

ht5_simple <- GroupHeatmap(ss2,assay = "RNA",
                    features = top_mk.hc$gene, feature_split = top_mk.hc$cluster, 
                    group.by = "merged_sub.anno_type",#split.by = "orig.ident", 
                    heatmap_palette = "YlOrRd",
                    #group_palcolor =as.character(npcols),# "#FFD9A4" "#FF8F6C" "#EEBDD1" "#9DE7D7" "#FFE671"
                    feature_split_palcolor=npcols,
                    #cell_split_palcolor = rev(c("#d53e4f","#fee08b","#abdda4","#3288bd")),
                    #group_palette ="viridis",
                    group_palcolor = list(c(npc)),column_title = "NP Cell clusters",
                      
                    #cell_annotation = c(nps),#,"Phase","foxp2","foxp4", "oxt"), 
                    #cell_annotation_palette = simcol,
                    #cell_annotation_palcolor = as.character(npcols),
                    #cell_annotation_params = list(height = unit(10, "mm")),
                    feature_annotation = c("TF","dGCprimed"),label_size = 1.5,#flip = T,
                    feature_annotation_palcolor = list(c("gold", "steelblue"),c("red")),
                    features_fontsize = c(6,25),dot_size = unit(3, "mm"),
                    add_dot = T, add_bg = TRUE, nlabel = 0, show_row_names = F
)

print(ht5_simple$plot)

pdf(paste0("./figures/cell_type/hcluster_top10nps_orig_mk_r4000_mergedsub.pdf"),
     width = 10,height = 25)
print(ht5_simple$plot)
dev.off()


#top5
top_mk.hc5= ss2.markers %>% 
  dplyr::filter(abs(avg_log2FC)>1) %>%
  dplyr::filter(p_val_adj <0.1) %>%
  dplyr::filter(pct.1 >0.3) %>%
  group_by(cluster) %>%
  top_n(5, wt= -p_val_adj )
  
ht5 <- GroupHeatmap(ss2,assay = "RNA",
                    features = top_mk.hc5$gene, feature_split = top_mk.hc5$cluster, 
                    group.by = "merged_sub.anno_type",#split.by = "orig.ident", 
                    heatmap_palette = "YlOrRd",
                    #group_palcolor ="",# "#FFD9A4" "#FF8F6C" "#EEBDD1" "#9DE7D7" "#FFE671"
                    feature_split_palcolor=as.list(npcols),
                    #cell_split_palcolor = rev(c("#d53e4f","#fee08b","#abdda4","#3288bd")),
                    group_palcolor = list(c(npc)),column_title = "NP Cell clusters",                    
                    #group_palette = "viridis",
                    #cell_annotation = c(nps),#,"Phase","foxp2","foxp4", "oxt"), 
                    #cell_annotation_palette = simcol,
                    #cell_annotation_params = list(height = unit(10, "mm")),
                    feature_annotation = c("TF","dGCprimed"),label_size = 1.5,#flip = T,
                    feature_annotation_palcolor = list(c("gold", "steelblue"),c("red")),
                    add_dot = T, add_bg = TRUE, nlabel = 0, show_row_names = TRUE
)

print(ht5$plot)

pdf(paste0("./figures/cell_type/hcluster_top5nps_orig_mk_r4000_mergedsub_wlabel.pdf"),
    width = 15,height = 20)
print(ht5$plot)
dev.off()

ht5_simple <- GroupHeatmap(ss2,assay = "RNA",
                           features = top_mk.hc5$gene, feature_split = top_mk.hc5$cluster, 
                           group.by = "merged_sub.anno_type",#split.by = "orig.ident", 
                           heatmap_palette = "YlOrRd",
                           #group_palcolor =as.character(npcols),# "#FFD9A4" "#FF8F6C" "#EEBDD1" "#9DE7D7" "#FFE671"
                           feature_split_palcolor=npcols,
                           #cell_split_palcolor = rev(c("#d53e4f","#fee08b","#abdda4","#3288bd")),
                           #group_palette ="viridis",
                           group_palcolor = list(c(npc)),column_title = "NP Cell clusters",
                           
                           #cell_annotation = c(nps),#,"Phase","foxp2","foxp4", "oxt"), 
                           #cell_annotation_palette = simcol,
                           #cell_annotation_palcolor = as.character(npcols),
                           #cell_annotation_params = list(height = unit(10, "mm")),
                           feature_annotation = c("TF","dGCprimed"),label_size = 1.5,#flip = T,
                           feature_annotation_palcolor = list(c("gold", "steelblue"),c("red")),
                           features_fontsize = c(6,25),dot_size = unit(3, "mm"),
                           add_dot = T, add_bg = TRUE, nlabel = 0, show_row_names = F
)

print(ht5_simple$plot)

pdf(paste0("./figures/cell_type/hcluster_top5nps_orig_mk_r4000_mergedsub.pdf"),
    width = 10,height = 20)
print(ht5_simple$plot)
dev.off()




saveRDS(
  ss2,
  file = paste0("./data/rds/step5_nps_cont_pre.var",r.variable,".rds")
)


######
#ATAC for npc
library(BSgenome.Drerio.UCSC.danRer11)
ss2=ss2.ori
DefaultAssay(ss2) <- "ATAC_macs3"

ss2 <- RegionStats(ss2, genome = BSgenome.Drerio.UCSC.danRer11)

# link peaks to genes
ss2 <- LinkPeaks(
  object = ss2,
  peak.assay = "ATAC_macs3",
  expression.assay = "RNA"
  
)


#atac-seq weight
library(Signac)


npcs_groups=unique(c(ss2$np50_clusters_mask[grep("_",ss2$np50_clusters_mask)]))
npc3=unique(ss2$merged_sub.anno_type[which(ss2$np50_clusters_mask %in%npcs_groups )])


library(presto)
library(GenomicRanges)
library(viridis)

DA_ct_npc3= wilcoxauc(ss2,
                      groups_use = npc3,
                      group_by = "merged_sub.anno_type", 
                      seurat_assay = "ATAC_macs3")

top_peaks_ct <- DA_ct_npc3 %>%
  dplyr::filter(abs(logFC) > log2(1.25) &
                  padj < 0.1
                #pct_in - pct_out > 13 &
                #auc > 0.55
  ) %>%
  group_by(group) 
ranges.show <- StringToGRanges(top_peaks_ct$feature)
ranges.show$color <- "gray"
subgroup2show=npc3

#DEGs
#markers

neuropep=c("oxt","avp","th","th2","sst1.1","galn","npvf","fshb","agrp","pmch","pomc",
           "hcrt","npffl","nmbb","nmba","edn1","edn2","edn3b","calca","prlh2","kiss1",
           "npy","pyya","vip","ccka","penka","penkb","nts","crhb","trh","tshba",
           "tshbb","gnrh3","gnrh2","lhb","cga")
top_nps= intersect(neuropep,ss2.markers$gene)

#nps
#####
dir.create("./figures/ATAC/")
dir.create("./figures/ATAC/nps")
dir.create("./figures/ATAC/markers")

for (i in top_nps) {
  gene_cordi=LookupGeneCoords(ss2,i)
  hits <- findOverlaps(ranges.show, gene_cordi)
  
  # Extract elements that fall within the range
  selected_elements <- ranges.show[queryHits(hits)]
  
  a=CoveragePlot(
    object = ss2,#split.by = "orig.ident",#peaks.group.by = "orig.ident",  
    group.by = "merged_sub.anno_type",  
    region = gene_cordi,
    features = i,
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
    pdf(paste0("./figures/ATAC/nps/ATAC_nps_",i,".pdf"),
         width = 8,height = 10)
    print(a & scale_fill_manual(values = magma(length(subgroup2show)+1,alpha = 0.5)))
    dev.off()
    aa=paste0(i," has enough peak")
    print(aa)
  }else{
    aa=paste0(i," not enough peak")
    print(aa)
  }
}
#markers
#####
top_mk_gene=top_mk.hc$gene

for (i in top_mk_gene) {
  gene_cordi=LookupGeneCoords(ss2,i)
  if(is.null(gene_cordi) == FALSE){
    hits <- findOverlaps(ranges.show, gene_cordi)
    
    # Extract elements that fall within the range
    selected_elements <- ranges.show[queryHits(hits)]
    
    a=CoveragePlot(
      object = ss2,#split.by = "orig.ident",#peaks.group.by = "orig.ident",  
      group.by = "merged_sub.anno_type",  
      region = gene_cordi,
      features = i,
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
      pdf(paste0("./figures/ATAC/markers/ATAC_maker_",i,".pdf"),
           width = 8,height = 10)
      print(a & scale_fill_manual(values = magma(length(subgroup2show)+1,alpha = 0.5)))
      dev.off()
      aa=paste0(i," has enough peak")
      print(aa)
    }else{
      aa=paste0(i," not enough peak")
      print(aa)
    }
  }
}






# save Seurat objects -----------------------------------------------------

#for nps
save(
  ss2.markers,
  top_mk.hc,
  npc,
  file = paste0("./data/rda/nps_markers_",r.variable,".rda")
)

saveRDS(
  ss2,
  file = paste0("./data/rds/step5_nps_var",r.variable,".rds")
)


