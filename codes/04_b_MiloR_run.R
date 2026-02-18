######
#milo for cell composition
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)

#r=4000
#from seurat
DefaultAssay(ss1)="SCT"

# Set identity (if needed)
Idents(ss1) <- "orig.ident"  # or use any metadata column

# Extract metadata
set.seed(2345)
meta_df <- ss1@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  group_by(orig.ident) %>%
  sample_n(size = 4500, replace = FALSE)  # Replace = TRUE if < 4000 cells in some groups

# Subset Seurat object with sampled cells
ss1_4500 <- subset(ss1, cells= meta_df$cell)

#format change into milo ss1_miloject
ss1_sce <- as.SingleCellExperiment(ss1_4500)
ss1_milo <- Milo(ss1_sce)

#build knn graph
ss1_milo <- buildGraph(ss1_milo, k = 25, d=100,reduced.dim = "HARMONY")
ss1_milo <- makeNhoods(ss1_milo, prop = 0.1, k = 25, d=100, refined = TRUE)
ss1_milo@nhoods
ss1_milo@nhoodIndex

cellIdx=sapply(ss1_milo@nhoodIndex, function(x){colnames(ss1_milo)[x[1]]})

#histo
plotNhoodSizeHist(ss1_milo)

#count cells in neighborhood
ss1_milo@colData$merged_sub.anno_type=gsub("\\/","",ss1_milo@colData$merged_sub.anno_type)
ss1_milo@colData$merged_sub.anno_type_ori=paste0(ss1_milo@colData$merged_sub.anno_type,".",ss1_milo@colData$orig.ident)

ss1_milo <- countCells(ss1_milo, meta.data = data.frame(colData(ss1_milo)), samples="merged_sub.anno_type_ori")
head(nhoodCounts(ss1_milo))

#####
#Differential abundance testing
#defining experimental design

a= str_split(ss1_milo@colData$orig.ident,"_")

geno=sapply(a, function(x){x[1]})
LD= sapply(a, function(x){x[2]})      

ss1_milo@colData$geno=geno
ss1_milo@colData$LD=LD

#####
### Using TMM normalisation
#Performing spatial FDR correction with k-distance weighting

  milo_design <- data.frame(colData(ss1_milo))[,c("merged_sub.anno_type_ori", "merged_sub.anno_type", "geno","LD","orig.ident")]
  milo_design <- distinct(milo_design)
  rownames(milo_design) <- milo_design$merged_sub.anno_type_ori
  ## Reorder rownames to match columns of nhoodCounts(milo)
  milo_design <- milo_design[colnames(nhoodCounts(ss1_milo)), , drop=FALSE]
  table(milo_design$orig.ident)
  
  rownames(milo_design) <- milo_design$merged_sub.anno_type_ori
  contrast.all <- c("orig.identcont_post - orig.identcont_pre","orig.identbPAC_post - orig.identbPAC_pre",
                    "orig.identbPAC_pre - orig.identcont_pre","orig.identbPAC_post - orig.identcont_post"
                    ) # the syntax is <VariableName><ConditionLevel> - <VariableName><ControlLevel>

# this is the edgeR code called by `testNhoods`
model <- model.matrix(~ 0 + orig.ident, data=milo_design)
mod.constrast <- makeContrasts(contrasts=contrast.all, levels=model)

mod.constrast

contrast1.res <- testNhoods(ss1_milo, design=~0+ orig.ident, design.df=milo_design, fdr.weighting="graph-overlap", model.contrasts = contrast.all)
head(contrast1.res)
table(contrast1.res$SpatialFDR < 0.1)

### compare specific point.
comparisons=c("contLD","bPACLD","preGeno","postGeno")
for (cp in 1:4) {
  name=paste0("contrast_",comparisons[cp])
  print(name)
  tn <- testNhoods(ss1_milo, design=~0+ orig.ident, design.df=milo_design, fdr.weighting="graph-overlap", 
                       model.contrasts = contrast.all[cp])
  head(tn)
  table(tn$SpatialFDR < 0.1)
  assign(name,tn)
}


#LogFC which compares nhood abundance 


a=plot(contrast1.res$logFC.orig.identcont_post...orig.identcont_pre, contrast_contLD$logFC,
     xlab="contPost vs. contPre LFC\nsingle contrast", ylab="contPost vs. contPre LFC\nmultiple contrast")

b=plot(contrast1.res$SpatialFDR, contrast_contLD$SpatialFDR,
     xlab="Spatial FDR\nsingle contrast", ylab="Spatial FDR\nmultiple contrast")

#by geno
model_geno <- model.matrix(~ 0 + geno, data=milo_design)
contrast.geno=c("genobPAC - genocont")
mod.constrast_geno <- makeContrasts(contrasts=contrast.geno, levels=model_geno)


contrast_geno <- testNhoods(ss1_milo, design=~0+ geno, design.df=milo_design, 
                            fdr.weighting="graph-overlap", model.contrasts = mod.constrast_geno)

head(contrast_geno)
table(contrast_geno$SpatialFDR < 0.1)


#by LD
model_LD <- model.matrix(~ 0 + LD, data=milo_design)
contrast.LD=c("LDpost - LDpre")
mod.constrast_LD <- makeContrasts(contrasts=contrast.LD, levels=model_LD)

contrast_LD <- testNhoods(ss1_milo, design=~0+ LD, design.df=milo_design, 
                            fdr.weighting="graph-overlap", model.contrasts = mod.constrast_LD)

head(contrast_LD)

table(contrast_LD$SpatialFDR < 0.1)


#####
##visualization- umap-like embedding

ss1_milo <- buildNhoodGraph(ss1_milo)

a=CellDimPlot(
  srt = ss1, group.by = "merged_sub.anno_type", seed = 0,
  reduction = "umap", theme_use = "theme_blank",
  label = T,label_insitu = T,label_repel = T
) 
b=plotNhoodGraphDA(ss1_milo, contrast1.res, alpha=0.1, min_size=5) #+plot_layout(guides="collect")
c=plotNhoodGraphDA(ss1_milo, contrast_contLD, alpha=0.1)+scale_fill_gradient2(low = "blue",mid = "lightyellow",high = "firebrick",breaks = c(-2, 0, 2)) #+plot_layout(guides="auto")
d=plotNhoodGraphDA(ss1_milo, contrast_bPACLD, alpha=0.1)+scale_fill_gradient2(low = "blue",mid = "lightyellow",high = "firebrick",breaks = c(-2, 0, 2)) #+plot_layout(guides="collect")
e=plotNhoodGraphDA(ss1_milo, contrast_preGeno, alpha=0.1)+scale_fill_gradient2(low = "blue",mid = "lightyellow",high = "firebrick",breaks = c(-2, 0, 2)) #+plot_layout(guides="collect")
f=plotNhoodGraphDA(ss1_milo, contrast_postGeno, alpha=0.1)+scale_fill_gradient2(low = "blue",mid = "lightyellow",high = "firebrick",breaks = c(-2, 0, 2)) #+plot_layout(guides="collect")
g=plotNhoodGraphDA(ss1_milo, contrast_geno, alpha=0.1, min_size=5) #+plot_layout(guides="collect")
h=plotNhoodGraphDA(ss1_milo, contrast_LD, alpha=0.1, min_size=5) #+plot_layout(guides="collect")

tiff(paste0("./figures/umap_RNA_all.tiff"),
     width = 30,height = 30,units = "cm", res = 300,compression = "lzw", bg = NA)
print((a))
dev.off()

tiff(paste0("./figures/milo_",r.variable,"_all.tiff"),
     width = 40,height = 20,units = "cm", res = 300,compression = "lzw", bg = NA)
print((a|b))
dev.off()


tiff(paste0("./figures/milo_LD_",r.variable,"_all.tiff"),
     width = 50,height = 40,units = "cm", res = 300,compression = "lzw", bg = NA)
print(((c|d)/(e|f)))
dev.off()

####more-visualization
comparisons=c("contLD","bPACLD","preGeno","postGeno")
contrast_name=paste0("contrast_",comparisons)

clinfo=ss1[["merged_sub.anno_type"]][cellIdx,]

cell_node=data.frame("cell_id"=cellIdx,
                     "cluster"=clinfo)
DA_signnodes=list()
sign_milo=list()

for (i in contrast_name) {
  tb=get(i)
  tb= cbind (tb, cell_node)
  #tb_ft=tb %>% dplyr::filter(SpatialFDR < 0.1) 
  tb$source=i
  sign_milo[[i]]=tb
  tb_ft2=tb %>% dplyr::filter(abs(logFC) > 1) 
  a=table(tb_ft2$cluster)
  DA_signnodes[[i]]=a  
}

DA_signnodes_tb=data.frame(DA_signnodes[[1]],DA_signnodes[[2]],
                           DA_signnodes[[3]],DA_signnodes[[4]])
DA_signnodes_tb=DA_signnodes_tb[,c(1:2,4,6,8)]
colnames(DA_signnodes_tb)=c("cluster",contrast_name)
write.csv(DA_signnodes_tb, paste0("./outputs/cell_type/milo_DA_signnodes_tb.csv"))

sign_milo_tb=rbind(sign_milo[[1]],sign_milo[[2]],sign_milo[[3]],sign_milo[[4]])

write.csv(sign_milo_tb, paste0("./outputs/cell_type/sign_milo_tb.csv"))

######
#visualization

# library
library(ggplot2)
library(ggExtra)

# data
head(sign_milo_tb)
sign_milo_tb_ft=dplyr::filter(sign_milo_tb,source %in% c("contrast_contLD","contrast_bPACLD"))
sign_milo_tb_ft$dir=sign_milo_tb_ft$logFC
sign_milo_tb_ft$dir[sign_milo_tb_ft$dir <0 ]= 0
sign_milo_tb_ft$dir[sign_milo_tb_ft$dir >0 ]= 1
sign_milo_tb_ft$dir=as.character(sign_milo_tb_ft$dir)
sign_milo_tb_ft$dir[sign_milo_tb_ft$dir == "0" ]= "down"
sign_milo_tb_ft$dir[sign_milo_tb_ft$dir == "1" ]= "up"
sign_milo_tb_ft$cluster=factor(sign_milo_tb_ft$cluster,levels = type_table_m.v2$numbered)

# Your significance threshold
threshold <- -log10(0.01)

sign_milo_tb_ft <- sign_milo_tb_ft %>%
  mutate(
    # Create a new column for coloring.
    # If a point is below the threshold, label it "Non-significant".
    # Otherwise, keep its original 'source' label.
    plot_color = ifelse(-log10(SpatialFDR) < threshold, "Non-significant", as.character(source))
  )

# classic plot :


p <- ggplot(sign_milo_tb_ft, aes(x=cluster, y= -log10(SpatialFDR), fill=plot_color,
                                 size=abs(logFC), shape = dir)) +
  
  # Points are unchanged, still with a black outline
  geom_point(position = position_dodge(width = 0.8), color = "gray", stroke = 0.5) +
  
  # The significance line
  geom_hline(aes(linetype = "gas guzzler", yintercept = -log10(0.01)), color = "#d95f02") +
  
  # Update the fill label to be more descriptive
  labs(fill="Comparison",linetype ="sFDR=0.01", size="abs(LogFC)", shape="Regulation") +
  
  theme(panel.background = element_blank(),strip.background = element_rect(colour="gray"),
        axis.text.x=element_text(color = "black", size=10, angle=90, vjust=.8, hjust=0.8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.major = element_line(colour = "gray")) +
  
  # --- KEY CHANGES ARE HERE ---
  # 1. Add your new "Non-significant" category with a grey color
  scale_fill_manual(values = c("contrast_contLD" = alpha("#440154FF",0.5), 
                               "contrast_bPACLD" = alpha("#FDE725FF",0.8), 
                               "Non-significant" = "grey80")) +
  
  # The shape scale remains the same
  scale_shape_manual(values = c("up" = 24, "down" = 25)) +
  # The legend fix from before still works perfectly
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5)))

p + coord_flip()


tiff(paste0("./figures/milo_LD_plot_v.tiff"),
     width = 13,height = 30,units = "cm", res = 300,compression = "lzw", bg = NA)
print(p + scale_x_discrete(limits = rev(levels(sign_milo_tb_ft$cluster)))+coord_flip())
dev.off()

tiff(paste0("./figures/milo_LD_plot_h.tiff"),
     width = 30,height = 10,units = "cm", res = 300,compression = "lzw", bg = NA)
print(p + theme(legend.text = element_text(size=8),legend.position = "top"))
dev.off()

#table
sign_milo_tb_ft2=sign_milo_tb %>% dplyr::filter(abs(logFC) > 1) %>% dplyr::filter(SpatialFDR < 0.01)
write.csv(sign_milo_tb_ft2,paste0("./outputs/sign_milo_tb_ft2.csv"))

clusters=names(table(sign_milo_tb_ft2$cluster))
count_cluster=list()
for (i in contrast_name) {
  tb= sign_milo_tb_ft2 %>% dplyr::filter(source == i)
  count=table(tb$cluster)
  count_cluster[[i]]=as.numeric(count)
}
count_cluster_tb=data.frame("cluster"=clusters,
                            "contrast_contLD"=count_cluster[[1]],
                            "contrast_bPACLD"=count_cluster[[2]],
                            "contrast_preGeno"=count_cluster[[3]],
                            "contrast_postGeno"=count_cluster[[4]])

write.csv(count_cluster_tb,paste0("./outputs/count_cluster_tb_sign_milo_tb_ft2.csv"))


sign_milo_tb <- sign_milo_tb %>%
  mutate(
    # Create a new column for coloring.
    # If a point is below the threshold, label it "Non-significant".
    # Otherwise, keep its original 'source' label.
    plot_color = ifelse(-log10(SpatialFDR) < -log10(0.01)|abs(logFC) < 1, "Non-significant", as.character(source))
  )
sign_milo_tb$cluster

##reorder colonm by hclust result

type_table_m.v2= read.csv(paste0("./outputs/cell_type/cell_type_table_m_modified.csv"), header = T, row.names = 1)
type_table_m.v3= type_table_m.v2[- c(grep(pattern = "UnD", type_table_m.v2$numbered)),]

load(file=paste0("./data/rda/hcluster.rda"))

taxa_order_tb=data.frame("label"=taxa_order[- c(grep("UnD",taxa_order))]) 
tmp=join_by("label"=="numbered_mk")
taxa_order_tb=left_join(taxa_order_tb,type_table_m.v3,by=tmp)

sign_milo_tb$cluster=factor(sign_milo_tb$cluster, levels =taxa_order_tb$numbered )


milotb <- ggplot() +
  # Background layer: non-significant points
  geom_point(data = subset(sign_milo_tb, plot_color == "Non-significant"),
             aes(x = cluster, y = logFC, size = -log10(SpatialFDR), 
                 shape = source, fill = plot_color),
             position = position_dodge(width = 0),
             color = "gray90") +  # optional gray border
  
  # Foreground layer: significant (colored) points
  geom_point(data = subset(sign_milo_tb, plot_color != "Non-significant"),
             aes(x = cluster, y = logFC, size = -log10(SpatialFDR), 
                 shape = source, fill = plot_color),
             position = position_dodge(width = 1),
             color = "gray30") +  # optional darker border
  
  # Horizontal lines
  geom_hline(yintercept = 1, linetype = "dotdash", color = alpha("#d95f02", 0.5)) +
  geom_hline(yintercept = -1, linetype = "dotdash", color = alpha("#d95f02", 0.5)) +
  
  # Theme and scales
  theme(panel.background = element_blank(),
        strip.background = element_rect(colour = "lightgray"),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, vjust = 0.8, hjust = 0.8),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.major = element_line(colour = "lightgray")) +
  
  scale_shape_manual(values = c(
    "contrast_bPACLD" = 21, 
    "contrast_contLD" = 22, 
    "contrast_postGeno" = 23,
    "contrast_preGeno" = 24
  )) +
  
  scale_fill_manual(values = c(
    "contrast_bPACLD" = alpha("#d53e4f", 0.7),
    "contrast_contLD" = alpha("#3288bd", 0.7),
    "contrast_postGeno" = alpha("#abdda4", 0.7),
    "contrast_preGeno" = alpha("#fee08b", 0.7),
    "Non-significant" = alpha("gray90", 0)
  )) +
  
  guides(
    shape = guide_legend(override.aes = list(
      fill = c(
        alpha("#d53e4f", 0.7),
        alpha("#3288bd", 0.7),
        alpha("#abdda4", 0.7),
        alpha("#fee08b", 0.7)
      ),
      shape = c(21, 22, 23, 24),
      size = 4
    )),
    
    fill = guide_legend(override.aes = list(
      shape = 21,
      size = 4
    ))
  ) +
  
  scale_x_discrete(limits = rev) + 
  coord_flip()


tiff(paste0("./figures/milo_LD_plot_v2.tiff"),
     width = 15,height = 30,units = "cm", res = 300,compression = "lzw", bg = NA)
print(milotb)
dev.off()

###save
save(cellIdx,sign_milo_tb,contrast1.res,contrast_contLD,contrast_bPACLD,contrast_preGeno,contrast_postGeno,contrast_geno,contrast_LD,
  file = paste0("./data/rda/ss1_milo.rda")
)

saveRDS(ss1_milo,paste0("./data/rds/ss1_milo.rds"))
