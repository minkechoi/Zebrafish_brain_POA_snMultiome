library(plyr)
library(dplyr)
library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(Matrix)
library(viridis)

fcha <- function(){ gsub("-","",Sys.Date()) }

dir <- "~/projects/smed_cisreg/"
setwd(dir)

outdir <- "~/projects/smed_cisreg/outputs/iterations_clustering/"

# Define custom normalization methods
do_pf <- function(mtx, sf = NULL) {
  pf <- Matrix(rowSums(mtx), ncol = 1) # sum of counts per cell
  if (is.null(sf)) {
    sf <- mean(pf) # calculate scaling factor
  }
  dsf <- Diagonal(x = as.numeric(sf/pf)) # diagonal sparse matrix with sf/pf values on the diagonal
  res <- dsf %*% mtx # matrix multiplication should return same dimensions as cts
  rownames(res) <- rownames(mtx)
  colnames(res) <- colnames(mtx)
  
  return(res)
}

# Read in all input expression matrices
print("loading dataset")
smed_cisreg_scrna <- readRDS("outputs/rda/20231114_smed_cisreg_scRNA_concatenated_libraries.RDS")

# Loop
npcs <- c(50,70,95,120)
hvgs <- c(15000,20000)
neighbours <- c(35,50,65)
resolutions <- c(1.5,2,2.5)

for (pc in npcs){
  for(hvg in hvgs){
    for(neighbour in neighbours){
      for(resolution in resolutions){
        
        print(paste0("Starting with ",hvg," highly variable genes, ",pc," PCs, ",neighbour," neighbours, resolution: ",resolution))
        
        scdata <- smed_cisreg_scrna
        
        # resampling variable genes
        print("finding variable genes")
        var.genes <- SelectIntegrationFeatures(SplitObject(scdata, split.by = "library"), nfeatures = hvg, verbose = TRUE, fvf.nfeatures = hvg, selection.method = "vst")
        VariableFeatures(scdata) <- var.genes
        
        scdata <- ScaleData(scdata, features = VariableFeatures(scdata))
        
        # Do PCA on data including only the variable genes.
        print("PCA")
        ndims <- ifelse(pc > 25,50,25)
        
        scdata <- RunPCA(
          scdata, 
          features = VariableFeatures(scdata), 
          npcs = pc, 
          ndims.print = 1:ndims, 
          nfeatures.print = 5
        )
        
        # Plotted overlapping
        dimplot <- DimPlot(scdata, reduction = "pca", dims = c(1, 2), group.by = "library", pt.size = 1, raster = FALSE)
        pdf(paste0(outdir,fcha(),"_smed_cisreg_scrna_",hvg,"_hvgs_",pc,"_pcs_",neighbour,"_neighbours_res_",resolution,"_pca.pdf"))
        print(dimplot)
        dev.off()
        
        # Harmony
        print("Running Harmony")
        scdata <- 
          RunHarmony(
            scdata, "orig.ident",
            dims.use = 1:pc, 
            theta = 3, 
            lambda = 3, 
            nclust = 40,
            max.iter.harmony = 20,
            plot_convergence = TRUE
          )
        
        # Cluster the cells
        print("Finding neighbours")
        scdata <- FindNeighbors(scdata, reduction = "harmony", dims = 1:pc, k.param = neighbour)
        print("Finding clusters")
        scdata <- FindClusters(scdata, resolution = resolution, algorithm = 1, random.seed = 75)
        
        # Create a UMAP visualization.
        print("Running UMAP")
        scdata <- RunUMAP(scdata, dims = 1:pc, reduction = "harmony", n.neighbors = neighbour, min.dist = 0.5, spread = 1, metric = "euclidean", seed.use = 1)  
        
        # Plot umap
        Umapbycluster <- DimPlot(scdata, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 1, raster = FALSE)
        pdf(paste0(outdir,fcha(), "_smed_cisreg_scrna_",hvg,"_hvgs_",pc,"_pcs_",neighbour,"_neighbours_res_",resolution,"_Umap_by_cluster.pdf"))
        print(Umapbycluster)
        dev.off()
        
        Umapbycluster <- DimPlot(scdata, reduction = "umap", group.by = "library", pt.size = 1, raster = FALSE)
        pdf(paste0(outdir,fcha(),"_smed_cisreg_scrna_",hvg,"_hvgs_",pc,"_pcs_",neighbour,"_neighbours_res_",resolution,"_Umap_by_lib.pdf"))
        print(Umapbycluster)
        dev.off()
        
        Umapbycluster <- DimPlot(scdata, reduction = "umap", split.by = "orig.ident", pt.size = 1, raster = FALSE)
        pdf(paste0(outdir,fcha(),"_smed_cisreg_scrna_",hvg,"_hvgs_",pc,"_pcs_",neighbour,"_neighbours_res_",resolution,"_Umap_by_ident.pdf"), width = 30, height = 8)
        print(Umapbycluster)
        dev.off()
        
        # print("Finding markers")
        # scdata.markers <- 
        #   FindAllMarkers(
        #     scdata,
        #     only.pos = TRUE, 
        #     return.thresh = 1,
        #     logfc.threshold = 0
        #   )
        # print("Done finding markers")
        
        # markers <- 
        #   scdata.markers[
        #     scdata.markers$p_val_adj < 0.01,
        #   ]
        # 
        # markers <- 
        #   scdata.markers[
        #     scdata.markers$avg_log2FC > 0.005,
        #   ]
        # 
        # markers$cluster <- as.factor(markers$cluster)
        # 
        # markers_top <- markers %>%
        #   group_by(cluster) %>%
        #   top_n(n = 30, wt = avg_log2FC)
        # 
        # write.csv(markers, file = paste0(outdir,fcha(),"_smed_cisreg_scrna_",hvg,"_hvgs_",pc,"_pcs_",neighbour,"_neighbours_res_",resolution,"_markers.csv"))
        # write.csv(markers_top, file = paste0(outdir,fcha(),"_smed_cisreg_scrna_",hvg,"_hvgs_",pc,"_pcs_",neighbour,"_neighbours_res_",resolution,"_markers_top.csv"))
        # write.csv(table(Idents(scdata)), file = paste0(outdir,fcha(),"_smed_cisreg_scrna_",hvg,"_hvgs_",pc,"_pcs_",neighbour,"_neighbours_res_",resolution,"_merged_numbers.csv"))
      
        # Save the dataset
        # print("Saving dataset")
        # saveRDS(
        #   scdata,
        #   file = paste0("outputs/rscript_outputs/",fcha(),"_smed_cisreg_scrna_",hvg,"_hvgs_",pc,"_pcs_",neighbour,"_neighbours_res_",resolution,"_dataset.rds")
        # )
        
        # Plot ncells per cluster
        print("plot num cells per cluster")
        # Create a data frame with the number of cells in each library and cluster
        df <- data.frame(table(scdata$orig.ident, scdata$seurat_clusters))
        
        # Rename the columns of the data frame
        colnames(df) <- c("Library", "Cluster", "Count")
        
        # Normalize the counts in each cluster
        df$logCount <- log(df$Count,10) / log(tapply(df$Count, df$Cluster, sum)[df$Cluster],10)
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
        
        # Show the grid of plots
        pdf(
          file = paste0(outdir,fcha(),"_smed_cisreg_scrna_",hvg,"_hvgs_",pc,"_pcs_",neighbour,"_neighbours_res_",resolution,"_numcells_cluster_library.pdf"),
          width = 10, height = 15
        )
        print(grid)
        dev.off()
        
        # Prepare plots
        print("various plots")
        #piwi
        piwi <- "h1SMcG0013999"
        scdata_piwi <-
          FeaturePlot(
            scdata,
            features = piwi,
            min.cutoff = 0, #this sets the limit for showing expression
            cols = c("lightgrey","#a100d2"),
            pt.size = 0.75,
            order = T,
            raster = FALSE
          )+ ggtitle("piwi")
        
        nanos <- "h1SMcG0003273"
        scdata_nanos <-
          FeaturePlot(
            scdata,
            features = nanos,
            min.cutoff = 0, #this sets the limit for showing expression
            cols = c("lightgrey","#a100d2"),
            pt.size = 0.75,
            order = T,
            raster = FALSE
          )+ ggtitle("nanos")
        
        estrella <- "h1SMcG0019080"
        scdata_estrella <-
          FeaturePlot(
            scdata,
            features = estrella,
            min.cutoff = 0, #this sets the limit for showing expression
            cols = c("lightgrey","#a100d2"),
            pt.size = 0.75,
            order = T, 
            raster = FALSE
          )+ ggtitle("estrella+ (glia)")
        
        # Big PDF
        pdf(
          file = paste0(outdir,fcha(),"_smed_cisreg_scrna_",hvg,"_hvgs_",pc,"_pcs_",neighbour,"_neighbours_res_",resolution,"_featplots.pdf"),
          width = 26,
          height = 10
        )
        print(
          plot_grid(
            scdata_piwi,
            scdata_nanos,
            scdata_estrella,
            ncol = 3
          )
        )
        dev.off()
        
        print(paste0("Done with ",neighbour," neighbours, resolution: ",resolution))
      } 
    }
  }
}