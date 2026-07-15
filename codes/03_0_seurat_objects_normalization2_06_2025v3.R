# =============================================================================
# 03_0_seurat_objects_normalization2_06_2025v3.R
# -----------------------------------------------------------------------------
# Purpose : Core normalisation, batch integration and joint RNA+ATAC clustering
#           of the QC-filtered multiome object. Steps: SCTransform + cell-cycle
#           regression (RNA), PCA, Harmony batch correction, RNA clustering/UMAP,
#           feature annotation (TF / GC-primed genes), ATAC LSI + Harmony, and
#           weighted-nearest-neighbour (WNN) joint clustering/UMAP.
# Inputs  : ./data/rds/ft_comb_seurat_macs.rds  (QC-filtered, MACS3 peaks)
#           ./data/regev_lab_cell_cycle_genes_asFish.txt (cell-cycle genes)
#           ./data/Table S4. List of Adult LD-DEGs.csv  (GC-primed gene list)
# Outputs : intermediate + final RDS under ./data/rds/ (final:
#           step3_norm_harmony_wnn_RNA_ft<r.variable>.rds), UMAP/QC figures.
# Note    : PC1 is deliberately dropped from Harmony/UMAP as it tracks depth.
# =============================================================================

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc()

library(Seurat)
library(Signac)
library(tidyverse)
library(BSgenome.Drerio.UCSC.danRer11)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
library(rtracklayer)
library(harmony)          # batch integration
library(colorspace)
library(RColorBrewer)

set.seed(1234)
options(future.globals.maxSize = 120000 * 1024^2)

#set working dir

r.variable=4000           # number of SCT variable features
vs="06_2025_v4"
setwd(paste0("./",vs))


#functions adopted and modified from Alberto Perez-Posada @apposada
source("./ext_code/r_code/functions/sourcefolder.R")
sourceFolder("./ext_code/r_code/functions/")

#other temp function: negated %in%
'%!in%' <- function(x,y)!('%in%'(x,y))


#####
#RNA assay

# normalization_and_weight ------------------------------------------------

#######load cell cycling data
# Import cell cycle genes and calculate cell cycle scores after normalization
# (Regev-lab S/G2M gene lists mapped to zebrafish orthologues)

fishCCgenes = readLines(con = "./data/regev_lab_cell_cycle_genes_asFish.txt")
s.genes = fishCCgenes[1:42]       # S-phase markers
g2m.genes = fishCCgenes[43:96]    # G2/M-phase markers

ft_comb_seurat=readRDS("./data/rds/ft_comb_seurat_macs.rds")


# finding variable features -----------------------------------------------

DefaultAssay(ft_comb_seurat)="RNA"

print(paste("Normalization and scaling, set var.feature =",r.variable))

# SCTransform and cell cycle regression
# (regress technical depth + mito%; keep top r.variable variable genes)
ft_comb_seurat = SCTransform(ft_comb_seurat, assay = "RNA",new.assay.name = "SCT",
                             vars.to.regress = c("nCount_RNA", "percent.mt"),
                             variable.features.n =r.variable, seed.use = 0,
                             verbose = FALSE)
var.genes=VariableFeatures(ft_comb_seurat)

# Score cell-cycle phase; CC.Difference lets us regress cycle while keeping
# proliferating-vs-quiescent signal (Seurat "alternative" workflow)
ft_comb_seurat = CellCycleScoring(ft_comb_seurat,
                                  s.features = s.genes,
                                  g2m.features = g2m.genes)

ft_comb_seurat$Phase = factor(ft_comb_seurat$Phase, levels = c('G1', 'S', 'G2M'))
ft_comb_seurat$CC.Difference = ft_comb_seurat$S.Score - ft_comb_seurat$G2M.Score



DefaultAssay(ft_comb_seurat)<-"RNA"

# Also produce a standard log-normalised + cell-cycle-regressed RNA assay
ft_comb_seurat=ft_comb_seurat %>%
  NormalizeData() %>%
  ScaleData(vars.to.regress = "CC.Difference",
            features = rownames(ft_comb_seurat))

DefaultAssay(ft_comb_seurat)<-"SCT"

BasicFeatures_presubset <-
  VlnPlot(ft_comb_seurat, features = c("nFeature_RNA", "nCount_RNA"))

print(BasicFeatures_presubset)


saveRDS(ft_comb_seurat, file=paste0("./data/rds/norm_RNA",r.variable,".rds"))


# PCA ---------------------------------------------------------------------

# Do PCA on data including only the variable genes.

DefaultAssay(ft_comb_seurat)="SCT"

ft_comb_seurat <-
  RunPCA(
    ft_comb_seurat, assay = "SCT",
    #features = var.genes,
    npcs = 150, set.seed = 0,
    ndims.print = 1:50,
    nfeatures.print = 5
  )

#list of the genes that define the top principal components:
Pca1_2 <-
  VizDimLoadings(
    ft_comb_seurat,
    dims = 1:4,
    reduction = "pca"
  )

print(Pca1_2)

# Plotted overlapping (dims 3 vs 2; PC1 avoided as it tracks library size)
dimplot <-
  DimPlot(
    ft_comb_seurat,
    reduction = "pca",
    dims = c(3, 2),
    group.by = "orig.ident",
    raster = FALSE
  )
print(dimplot)


# Harmony -----------------------------------------------------------------
# batch effect correction of all the libraries.
# (dims.use starts at 2 to exclude PC1 / depth component)

ft_comb_seurat <-
  RunHarmony(
    ft_comb_seurat, "orig.ident",
    dims.use = 2:150,
    theta = 2,
    lambda = 2,
    nclust = 45,
    max.iter.harmony = 20,
    plot_convergence = TRUE
  )

harmony_embeddings <- Embeddings(ft_comb_seurat, "harmony")

#cluster the cells using all of the PCs calculated above.
#This will embed cells in a knn-graph structure that can be helpful to identify data communities.

ft_comb_seurat <-
  FindNeighbors(
    ft_comb_seurat,assay = "SCT",
    reduction = "harmony",
    dims = 1:149,
    k.param = 35
  )


ft_comb_seurat <-
  FindClusters(
    ft_comb_seurat,
    cluster.name = "RNA.clusters",
    resolution = 3,
    algorithm = 4, #Leiden
    random.seed = 453
  )


# QC: RNA depth/complexity per cluster, split by sample
a=VlnPlot(ft_comb_seurat,group.by = "RNA.clusters",split.by = "orig.ident",
          cols = rev(c("blue","lightgreen","orange","red")),
          features = c("nCount_RNA","nFeature_RNA"),ncol = 1,
          pt.size = 0)
a


tiff(paste0("./figures/cell_distribution_by_cluster",r.variable,".tiff"),width = 80,height = 30,units = "cm",
     res = 300,compression = "lzw")
print(a)
dev.off()

####


#Create a few UMAP visualizations
#PCA1 dim hasn't been used to minimize library size effects


ft_comb_seurat <-
  RunUMAP(
    ft_comb_seurat,
    dims = 1:120,
    reduction = "harmony",
    n.neighbors = 35,
    min.dist = 0.3,
    spread = 1,
    metric = "euclidean", #"cosine", #"euclidean",
    seed.use = 543,
    n.components = 2,alpha = 1, gamma = 1.0
  )

#assign color for clusters

library(magrittr)
library(RColorBrewer)

# Assemble a large qualitative palette to cover many clusters
cl_colors <-
  c(divergingx_hcl(8,"ArmyRose"),
    divergingx_hcl(11,"RdYlBu"),
    divergingx_hcl(7,"Zissou 1"),
    divergingx_hcl(11,"Spectral"),
    divergingx_hcl(8,"Fall"),
    sequential_hcl(8,"Hawaii"),
    divergingx_hcl(8,"Cividis")
  )

num_clusters <- length(unique(ft_comb_seurat$RNA.clusters))

set.seed(368)
cols <- sample(unname(cl_colors),num_clusters)

#the visualisation of the UMAP after library integration using Harmony,
#with cells from the different libraries in different colors.

Umap_group_library <-
  DimPlot(
    ft_comb_seurat,
    reduction = "umap",
    group.by = "orig.ident",
    cols = brewer.pal(4,"Spectral"),
    pt.size = 1,
    raster = FALSE
  )

print(Umap_group_library)

# UMAP with the clusters
Umap_cluster <-
  DimPlot(
    ft_comb_seurat,
    cols = cols,
    reduction = "umap",
    group.by = "RNA.clusters",
    label = TRUE,
    raster = FALSE
  )
print(Umap_cluster)

#split origin
Umap_split_library <-
  DimPlot(
    ft_comb_seurat,
    cols = cols,
    reduction = "umap",
    split.by = "orig.ident",
    ncol = 2,
    raster = FALSE
  )
print(Umap_split_library)


#save plots
dir.create(paste0("./figures/umaps"))

tiff(paste0("./figures/umaps/umap_cluster_var.genes_",r.variable,".tiff"),
     width = 25,height = 50,units = "cm", res = 300,compression = "lzw")
print(Umap_group_library/Umap_cluster/Umap_split_library)
dev.off()


# Marker genes of interest (POA neuropeptides + stress/monoamine markers)
goi=c("oxt","avp","crhb","sst1.1",
      "foxp2","nr3c2","fosab","fkbp5",
      "galn","th","th2","trh")

p4 <- FeaturePlot(ft_comb_seurat, #split.by = "orig.ident",
                  features = goi,
                  reduction="umap",
                  order = T,pt.size = 0.7,
                  min.cutoff = 0.8,cols = c("grey80", "firebrick"),
                  ncol=4) & NoAxes() & NoLegend()

dir.create(paste0("./figures/umaps/featureplots"))

tiff(paste0("./figures/umaps/featureplots/umap_featureplots_",r.variable,".tiff"),
     width = 40,height = 30,units = "cm", res = 300,compression = "lzw")
print(p4)
dev.off()


# additional meta-data ----------------------------------------------------
# Flag genes that are "GC-primed" (from adult LD-DEG table) and annotate TFs,
# stored on the RNA and SCT feature metadata for later use.

#primed gene_load
adult_LD_DEG = vroom::vroom("./data/Table S4. List of Adult LD-DEGs.csv")
primed_gene_table=adult_LD_DEG %>% dplyr::filter(GC_primed == "yes")
primed_genes=unique(primed_gene_table$zfin_id_symbol)
library(SCP)
primed.g=row.names(ft_comb_seurat[["RNA"]]@meta.features)
primed.g[which(primed.g %in%primed_genes)]="primed"
primed.g[which(primed.g != "primed")]=NA
ft_comb_seurat[["RNA"]]@meta.features["dGCprimed"]=primed.g
ft_comb_seurat <- AnnotateFeatures(ft_comb_seurat,assays = "RNA",
                                   species = "Danio_rerio",Ensembl_version=111, db = c("TF"))

primed.g=row.names(ft_comb_seurat[["SCT"]]@meta.features)
primed.g[which(primed.g %in%primed_genes)]="primed"
primed.g[which(primed.g != "primed")]=NA
ft_comb_seurat[["SCT"]]@meta.features["dGCprimed"]=primed.g
ft_comb_seurat <- AnnotateFeatures(ft_comb_seurat,assays = "SCT",
                                   species = "Danio_rerio",Ensembl_version=111, db = c("TF"))


# wnn_RNA-ATAC ------------------------------------------------------------
# Joint RNA+ATAC analysis: process ATAC (LSI), Harmony-correct it, then combine
# both modalities with weighted-nearest-neighbours (WNN).

# ATAC peak assay: TF-IDF normalise, select features, LSI (SVD)
DefaultAssay(ft_comb_seurat) = "ATAC"

ft_comb_seurat = RunTFIDF(ft_comb_seurat,assay = "ATAC")
ft_comb_seurat = FindTopFeatures(ft_comb_seurat, min.cutoff = 'q0')
ft_comb_seurat = RunSVD(ft_comb_seurat,n = 150)

# Batch-correct the ATAC LSI embedding as well
ft_comb_seurat <-
  RunHarmony(
    ft_comb_seurat, "orig.ident",
    assay.use="ATAC",
    dims.use = 2:150,
    theta = 2,
    lambda = 2,
    nclust = 45,
    max.iter.harmony = 20,
    reduction.save="harmony.atac"
  )

harmony_embeddings.atac <- Embeddings(ft_comb_seurat, "harmony.atac")


# Weighted nearest neighbours: learn per-cell RNA vs ATAC weighting
ft_comb_seurat <- FindMultiModalNeighbors(
  object = ft_comb_seurat,
  reduction.list = list("harmony", "harmony.atac"),
  dims.list = list(1:50, 1:50),
  #modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# Joint (WNN) UMAP embedding
ft_comb_seurat <- RunUMAP(ft_comb_seurat,
                          nn.name = "weighted.nn",
                          n.neighbors = 25,seed.use = 753,
                          spread = 1,
                          metric = "cosine", #euclidean",#"cosine",
                          min.dist = 0.5,alpha = 1, gamma = 1.0,
                          reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

# Joint (WNN) Leiden clustering -> primary cluster labels used downstream
ft_comb_seurat <-
  FindClusters(
    ft_comb_seurat,
    graph.name = "wsnn",
    cluster.name = "wnn_cluster",
    resolution = 2,
    algorithm = 4, #Leiden
    random.seed = 453
  )
ft_comb_seurat$seurat_clusters=ft_comb_seurat$wnn_cluster

num_clusters <- length(unique(ft_comb_seurat$wnn_cluster))

set.seed(368)
cols <- sample(unname(cl_colors),num_clusters)

#plot: WNN clusters on the WNN UMAP vs the RNA UMAP
a=DimPlot(ft_comb_seurat, reduction = "wnn.umap", group.by = "wnn_cluster",
        cols = cols,
        label = TRUE, label.size = 2.5, repel = TRUE)
a
b=DimPlot(ft_comb_seurat, reduction = "umap", group.by = "wnn_cluster",
          cols = cols,
          label = TRUE, label.size = 2.5, repel = TRUE)


tiff(paste0("./figures/umaps/umap_wnn_",r.variable,".tiff"),
     width = 50,height = 20,units = "cm", res = 300,compression = "lzw")
print(a|b)
dev.off()

#plot: RNA depth per WNN cluster

c=VlnPlot(ft_comb_seurat,group.by = "wnn_cluster",split.by = "orig.ident",
          cols = rev(c("blue","lightgreen","orange","red")),
          features = c("nCount_RNA","nFeature_RNA"),ncol = 1,
          pt.size = 0)
c


tiff(paste0("./figures/cell_distribution_by_wnncluster",r.variable,".tiff"),width = 80,height = 30,units = "cm",
     res = 300,compression = "lzw")
print(c)
dev.off()

# ATAC depth per WNN cluster
d=VlnPlot(ft_comb_seurat,group.by = "wnn_cluster",split.by = "orig.ident",
          cols = rev(c("blue","lightgreen","orange","red")),
          features = c("nCount_ATAC","nFeature_ATAC"),ncol = 1,
          pt.size = 0)
d


tiff(paste0("./figures/cell_distribution_by_wnncluster_ATAC_",r.variable,".tiff"),width = 80,height = 30,units = "cm",
     res = 300,compression = "lzw")
print(d)
dev.off()

# Marker genes on the WNN UMAP
DefaultAssay(ft_comb_seurat) = "SCT"
p4 <- FeaturePlot(ft_comb_seurat, #split.by = "orig.ident",
                  features = goi,
                  reduction="wnn.umap",
                  order = T,pt.size = 1,
                  min.cutoff = 0.8,cols = c("grey80", "firebrick"),
                  ncol=4) & NoAxes() & NoLegend()

tiff(paste0("./figures/umaps/featureplots/umap_wnn_featureplots_var.genes_",r.variable,".tiff"),
     width = 40,height = 30,units = "cm", res = 300,compression = "lzw")
print(p4)
dev.off()



#umi_wnn umap : total RNA counts on the WNN UMAP
p5=FeatureDimPlot(
  srt = ft_comb_seurat, features = "nCount_RNA",
  assay = "SCT",
  #label_repel = T,label_repulsion = 50,pt.size = 0.7,
  palette = "viridis",
  seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, #ncol = 4,add_density = T,
  reduction = "wnn.umap", theme_use = "theme_blank"
)& NoLegend()

dir.create(paste0("./figures/cell_type"))
tiff(paste0("./figures/cell_type/var_",r.variable,"umI_wnn.umap.tiff"),
     width = 10,height = 10,units = "cm", res = 300,compression = "lzw",bg = NA)
print(p5)
dev.off()

dir.create(paste0("./data/rds"))
dir.create(paste0("./data/rda"))

saveRDS(ft_comb_seurat,file=paste0("./data/rds/step3_norm_harmony_wnn_RNA",r.variable,".rds"))

####remove_ low cell cluster
#cell number < 100  (keep the 49 well-populated WNN clusters)
ft_comb_seurat= subset(ft_comb_seurat,wnn_cluster %in% c(1:49))
saveRDS(ft_comb_seurat,file=paste0("./data/rds/step3_norm_harmony_wnn_RNA_ft",r.variable,".rds"))




