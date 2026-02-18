library(Seurat)
library(Signac)
library(AUCell)
library(SCP)
library(patchwork)

#set the data and annotations 

set.seed(1234)
options(future.globals.maxSize = 120000 * 1024^2)

#set working dir

r.variable=4000
vs="06_2025_v4"
#setwd(paste0("./",vs))

#other temp function
'%!in%' <- function(x,y)!('%in%'(x,y))


#####
#data loading
#ss1= readRDS(paste0("./data/rds/norm_step4_var",r.variable,".rds")) 

Defaultassay(ss)="RNA"
#gene sets
ergs = c("fosaa","fosab","fosb","fosl2","itm2cb","egr1","egr2a","egr2b","egr3","ier2a")
GRs = c("nr3c1","nr3c2")
fkbp5 = c("fkbp5")
allsig=unique(c(ergs,GRs,fkbp5))

signatures <- list(
  ergs = c("fosaa","fosab","fosb","fosl2","itm2cb","egr1","egr2a","egr2b","egr3","ier2a"),
  GRs = c("nr3c1","nr3c2"),
  fkbp5 = c("fkbp5"),
  allsig=unique(c(ergs,GRs,fkbp5))
)

#AUC calculation

#####AUC_score

AUCellfunc <- function(exprMatrix,genes){
  
  out <- tryCatch(
    {
      cells_rankings <- AUCell_buildRankings(exprMatrix,plotStats=FALSE)
      geneSets <- list(geneSet1=genes)
      cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
      cellsAUCretrieve = getAUC(cells_AUC)
      cells_AUCellScores = t(cellsAUCretrieve)
      cells_AUCellScores = data.frame(rownames(cells_AUCellScores),cells_AUCellScores)
      colnames(cells_AUCellScores) = c('SampleID','AUCell')
      row.names(cells_AUCellScores)= NULL
      return(cells_AUCellScores)
      
    },
    
    error = function(cond) {
      cells_AUCellScores = data.frame(SampleID = "NA", AUCell = "NA")
      return(cells_AUCellScores)
    }
  )
  return(out)
}


exp.mtx=ss1@assays$RNA$scale.data
cells_rankings <- AUCell_buildRankings(exp.mtx, plotStats=TRUE)

auc_ergs=AUCellfunc(ss1@assays$RNA$scale.data, signatures[[1]])
auc_GRs=AUCellfunc(ss1@assays$RNA$scale.data, signatures[[2]])
auc_fkbp5=AUCellfunc(ss1@assays$RNA$scale.data, signatures[[3]])
auc_all=AUCellfunc(ss1@assays$RNA$scale.data, signatures[[4]])

ss1$AUC_ergs=auc_ergs$AUCell
ss1$AUC_GRs=auc_GRs$AUCell
ss1$AUC_fkbp5=auc_fkbp5$AUCell
ss1$AUC_allsig=auc_all$AUCell

#manual
ergs_geneSets <- list(geneSet1=signatures[[1]])
GRgenes_geneSets <- list(geneSet1=signatures[[2]])
fkbp5_geneSets <- list(geneSet1=signatures[[3]])
all_geneSets <- list(geneSet1=signatures[[4]])

cells_AUC1 <- AUCell_calcAUC(ergs_geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_AUC2 <- AUCell_calcAUC(GRgenes_geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_AUC3 <- AUCell_calcAUC(fkbp5_geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_AUC4 <- AUCell_calcAUC(all_geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)


set.seed(333)

cells_assignment1 <- AUCell_exploreThresholds(cells_AUC1, plotHist=F, assign=TRUE) 
cells_assignment2 <- AUCell_exploreThresholds(cells_AUC2, plotHist=F, assign=TRUE) 
cells_assignment3 <- AUCell_exploreThresholds(cells_AUC3, plotHist=F, assign=TRUE) 
cells_assignment4 <- AUCell_exploreThresholds(cells_AUC4, plotHist=F, assign=TRUE) 

#numgene_cell
nGenesPerCell <- apply(ss1@assays$RNA$scale.data, 2, function(x) sum(x>0))

colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
cellColorNgenes <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(nGenesPerCell,breaks=10, right=FALSE,include.lowest=TRUE))], names(nGenesPerCell))
ss1$nGenesPerCell=nGenesPerCell
ss1$nGenesPerCell_col=cellColorNgenes

##
ss1$orig.ident=factor(ss1$orig.ident, levels = c("cont_pre",  "cont_post", "bPAC_pre",  "bPAC_post"))

'
a=AUCell_plotHist(cellsAUC = cells_AUC1[1,],ylim = c(0,500),
                  aucThr = cells_assignment1$geneSet1$aucThr$selected)
b=AUCell_plotHist(cellsAUC = cells_AUC2[1,],ylim = c(0,500),
                  aucThr = cells_assignment1$geneSet1$aucThr$selected)
c=AUCell_plotHist(cellsAUC = cells_AUC3[1,],ylim = c(0,500),
                  aucThr = cells_assignment1$geneSet1$aucThr$selected)
'
p01=FeatureDimPlot(
  srt = ss1, features = paste0("AUC_",names(signatures))[1],
  lower_cutoff = cells_assignment1$geneSet1$aucThr$selected,
  assay = "RNA",
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor =c("lightgray","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, 
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",title = "IEGs",
  theme_args=theme(strip.text = element_text(size = 20,face = c("bold.italic")) )
) 
p02=FeatureDimPlot(
  srt = ss1, features = paste0("AUC_",names(signatures))[2],
  lower_cutoff = cells_assignment2$geneSet1$aucThr$selected,
  assay = "RNA",
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor = c("lightgray","#ffffcc","#d9f0a3","#addd8e","#238443","#005a32"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE,
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",title = "GRs",
  theme_args=theme(strip.text =  element_text(size = 20,face = c("bold.italic")) )
  
)


p03=FeatureDimPlot(
  srt = ss1, features = paste0("AUC_",names(signatures))[3],
  lower_cutoff = cells_assignment3$geneSet1$aucThr$selected,
  assay = "RNA",
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor = c("lightgray","#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, 
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",title = "fkbp5",
  theme_args=theme(strip.text = element_text(size = 20,face = c("bold.italic")) )
  
)



p1=FeatureDimPlot(
  srt = ss1, features = paste0("AUC_",names(signatures))[1],
  lower_cutoff = cells_assignment1$geneSet1$aucThr$selected,
  assay = "RNA",split.by = "orig.ident", 
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor =c("lightgray","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, ncol = 4,
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",
  theme_args=theme(strip.text = element_text(size = 0))
  
)


p2=FeatureDimPlot(
  srt = ss1, features = paste0("AUC_",names(signatures))[2],
  lower_cutoff = cells_assignment2$geneSet1$aucThr$selected,
  assay = "RNA",split.by = "orig.ident", 
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor = c("lightgray","#ffffcc","#d9f0a3","#addd8e","#238443","#005a32"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, ncol = 4,
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",
  theme_args=theme(strip.text = element_text(size = 0))
  
)


p3=FeatureDimPlot(
  srt = ss1, features = paste0("AUC_",names(signatures))[3],
  lower_cutoff = cells_assignment3$geneSet1$aucThr$selected,
  assay = "RNA",split.by = "orig.ident", 
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor = c("lightgray","#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, ncol = 4,
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",
  theme_args=theme(strip.text = element_text(size = 0))
  
)


#show_palettes()
p04=FeatureDimPlot(
  srt = ss1, features = "nGenesPerCell",
  assay = "RNA",
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palette = "cividis" ,
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, 
  #add_density = T,
  reduction = umap, theme_use = "theme_blank"
)
p4=FeatureDimPlot(
  srt = ss1, features = "nGenesPerCell",
  assay = "RNA",split.by = "orig.ident", 
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palette = "cividis" ,
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, ncol = 4,
  #add_density = T,
  reduction = umap, theme_use = "theme_blank"
)


#######

type_table_m.v3= type_table_m.v2[- c(grep(pattern = "UnD", type_table_m.v2$numbered)),]

simcol=type_table_m.v3$broader_col
names(simcol)=type_table_m.v3$broader
p.cl.broad=DimPlot(ss1,group.by = "merged_sub.anno", seed = 0,
           reduction = umap,cols = simcol,alpha = 0.8,#pt.size = 0.7,
           label = F, label.size =3,stroke.size = 0.5,
           order=T
)+theme_blank()+ ggtitle(NULL) 

p.cl.broad.sp=DimPlot(ss1,group.by = "merged_sub.anno", seed = 0,split.by = "orig.ident",
           reduction = umap,cols = simcol,alpha = 0.8,pt.size = 0.7,
           label = F, label.size =3,stroke.size = 0.5,
           order=T
)+theme_blank()+ ggtitle(NULL)  & NoLegend()


cycol=type_table_m.v3$clcol
names(cycol)=type_table_m.v3$cl_numb
p.cl=DimPlot(ss1,group.by = "merged_sub_numb", seed = 0,
            reduction = umap,cols = cycol,alpha = 0.8,#pt.size = 0.7,
            label = T, label.size =3,stroke.size = 0.5,
            order=T
)+theme_blank()+ ggtitle(NULL) & NoLegend()
p.cl.sp=DimPlot(ss1,group.by = "merged_sub_numb", seed = 0,split.by = "orig.ident",
            reduction = umap,cols = cycol,alpha = 0.8,pt.size = 0.7,
            label = F, label.size =3,stroke.size = 0.5,
            order=T
)+theme_blank()+ ggtitle(NULL) & NoLegend()



tiff(paste0("./figures/AUC_score_",r.variable,".tiff"),
     width = 50,height = 60,units = "cm", res = 300,compression = "lzw")
print(((p.cl.broad|p.cl.broad.sp)+plot_layout(widths = c(1, 4)))/((p.cl|p.cl.sp)+plot_layout(widths = c(1, 4)))/(p01|p1)/(p02|p2)/(p03|p3)/(p04|p4))
dev.off()

p04=FeatureDimPlot(
  srt = ss1, features = paste0("AUC_",names(signatures))[4],
  lower_cutoff = cells_assignment4$geneSet1$aucThr$selected,
  assay = "RNA",
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor = c("lightgray","#ffffb2","#e1bee7","#ba68c8","#aa00ff","#6a1b9a"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, 
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",title = "combined",
  theme_args=theme(strip.text = element_text(size = 20,face = c("bold.italic")) )
  
)



p4=FeatureDimPlot(
  srt = ss1, features = paste0("AUC_",names(signatures))[4],
  lower_cutoff = cells_assignment4$geneSet1$aucThr$selected,
  assay = "RNA",split.by = "orig.ident", 
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor =c("lightgray","#ffffb2","#e1bee7","#ba68c8","#aa00ff","#6a1b9a"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, ncol = 4,
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",
  theme_args=theme(strip.text = element_text(size = 0))
  
)

tiff(paste0("./figures/AUC_all_score_",r.variable,".tiff"),
     width = 50,height = 10,units = "cm", res = 300,compression = "lzw")
print((p04|p4))
dev.off()

####UCell
library(UCell)

#score calculation
#Ucell
ss1 <- AddModuleScore_UCell(ss1, 
                            features=signatures, name="Ucell_score")

ss1 <- SmoothKNN(ss1,
                 signature.names = paste0(names(signatures),"Ucell_score"),
                 reduction="harmony",k=20, suffix = "_smooth_ucell")


#visualization
#highlight cells from 35 and 45
c3545=names(ss1$merged_sub_numb[grepl("35|45",ss1$merged_sub_numb)])
p01=FeatureDimPlot(
  srt = ss1, features = paste0(names(signatures),"Ucell_score_smooth_ucell")[1],
  assay = "RNA",cells.highlight = c3545,
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor =c("lightgray","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, 
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",title = "IEGs",
  theme_args=theme(strip.text =  element_text(size = 0,face = c("bold.italic")) )
)
p02=FeatureDimPlot(
  srt = ss1, features = paste0(names(signatures),"Ucell_score_smooth_ucell")[2],
  assay = "RNA",cells.highlight = c3545,
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor = c("lightgray","#ffffcc","#d9f0a3","#addd8e","#238443","#005a32"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, 
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",title = "GRs",
  theme_args=theme(strip.text =  element_text(size = 0,face = c("bold.italic")) )
)


p03=FeatureDimPlot(
  srt = ss1, features = paste0(names(signatures),"Ucell_score_smooth_ucell")[3],
  assay = "RNA",cells.highlight = c3545,
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor = c("lightgray","#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, 
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",title = "fkbp5",
  theme_args=theme(strip.text =  element_text(size = 0,face = c("bold.italic")) )
)
p04=FeatureDimPlot(
  srt = ss1, features = paste0(names(signatures),"Ucell_score_smooth_ucell")[4],
  assay = "RNA",cells.highlight = c3545,
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor = c("lightgray","#ffffb2","#e1bee7","#ba68c8","#aa00ff","#6a1b9a"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, 
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",title = "combined",
  theme_args=theme(strip.text =  element_text(size = 0,face = c("bold.italic")) )
)




p1=FeatureDimPlot(
  srt = ss1, features = paste0(names(signatures),"Ucell_score_smooth_ucell")[1],
  assay = "RNA",split.by = "orig.ident", 
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor =c("lightgray","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, ncol = 4,
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",
  theme_args=theme(strip.text =  element_text(size = 0) )
)


p2=FeatureDimPlot(
  srt = ss1, features = paste0(names(signatures),"Ucell_score_smooth_ucell")[2],
  assay = "RNA",split.by = "orig.ident", 
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor = c("lightgray","#ffffcc","#d9f0a3","#addd8e","#238443","#005a32"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, ncol = 4,
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",
  theme_args=theme(strip.text =  element_text(size = 0) )
)


p3=FeatureDimPlot(
  srt = ss1, features = paste0(names(signatures),"Ucell_score_smooth_ucell")[3],
  assay = "RNA",split.by = "orig.ident", 
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor = c("lightgray","#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, ncol = 4,
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",
  theme_args=theme(strip.text =  element_text(size = 0) )
)


p4=FeatureDimPlot(
  srt = ss1, features = paste0(names(signatures),"Ucell_score_smooth_ucell")[4],
  assay = "RNA",split.by = "orig.ident", 
  #label_repel = T,label_repulsion = 50,
  pt.size = 0.7,palcolor =c("lightgray","#ffffb2","#e1bee7","#ba68c8","#aa00ff","#6a1b9a"),
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, ncol = 4,
  #add_density = T,
  reduction = umap, theme_use = "theme_blank",
  theme_args=theme(strip.text =  element_text(size = 0) )
)


#numgene_cell

pdf(paste0("./figures/UCell_score_",r.variable,".pdf"),
     width = 25,height = 30)
print(((p.cl.broad|p.cl.broad.sp)+plot_layout(widths = c(1, 4)))/((p.cl|p.cl.sp)+plot_layout(widths = c(1, 4)))/(p01|p1)/(p02|p2)/(p03|p3)/(p04|p4))
dev.off()

pdf(paste0("./figures/fig_2A-D_UCell_score_",r.variable,".pdf"),
    width = 27,height = 6)
print((p01|p02|p03|p04))
dev.off()

### score table 

Ucell_scores=ss1[[c("merged_sub.anno_type",
                   "ergsUcell_score_smooth_ucell","GRsUcell_score_smooth_ucell","fkbp5Ucell_score_smooth_ucell",
                   "AUC_ergs","AUC_GRs","AUC_fkbp5",
                   "orig.ident")
                  ]]

Ucell_scores$cl.orig=paste0(ss1$merged_sub.anno_type,":",ss1$orig.ident)


df= Ucell_scores[,c(1:7)]%>% group_by(merged_sub.anno_type) %>% summarise_all(mean)
write.csv(df,paste0("./outputs/AUCUCell_scores_by_cluster.csv"))

df= Ucell_scores[,c(9,2:7)]%>% group_by(cl.orig) %>% summarise_all(mean)

write.csv(df,paste0("./outputs/AUCUCell_scores_by_clusterand origin.csv"))

