# =============================================================================
# scType.R
# -----------------------------------------------------------------------------
# Purpose : Marker-based, reference-free cell-type annotation of clusters using
#           ScType (Ianevski et al., 2022). Scores each nucleus against curated
#           brain marker gene sets and assigns the top-scoring type per cluster.
# Context : Applied to the Seurat object `ss` (expects a per-cell cluster label
#           in `ss@meta.data$merged_sub` and a scaled RNA matrix).
# Output  : `sctype_scores` — one predicted type per cluster (low-confidence
#           clusters relabelled "Unknown").
# Note    : ScptType helper functions and the marker DB are sourced from the
#           upstream ScType GitHub repository (external dependency).
# =============================================================================

# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
# load gene set preparation function (ScType helper)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function (ScType scoring)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# DB file: ScType's curated marker database; select the Brain tissue markers
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus

# prepare gene sets: positive/negative marker lists for each brain cell type
gs_list <- gene_sets_prepare(db_, tissue)

# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(ss[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix (slot access differs between v4 and v5)
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(ss[["RNA"]]$scale.data) else as.matrix(ss[["RNA"]]@scale.data)



# run RNAype: per-cell enrichment scores for each candidate cell type
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. For raw (unscaled) count matrix set scaled = FALSE
# When using Seurat, we use "RNA" slot with 'scale.data' by default. Please change "RNA" to "RNA" for RNAransform-normalized data,
# or to "integrated" for joint dataset analysis. To apply RNAype with unscaled data, use e.g. ss[["RNA"]]$counts or ss[["RNA"]]@counts, with scaled set to FALSE.

# merge by cluster: sum per-cell scores within each cluster, keep top 10 types
cL_resutls <- do.call("rbind", lapply(unique(ss@meta.data$merged_sub), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(ss@meta.data[ss@meta.data$merged_sub==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(ss@meta.data$merged_sub==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  # winning type per cluster

# set low-confident (low RNAype score) clusters to "unknown"
# (confidence heuristic: top score must exceed ncells/4)
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])

