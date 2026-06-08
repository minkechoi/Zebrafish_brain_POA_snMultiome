# must be called from Rmarkdown
fragments = c(500,1000)
pcts = c(10,15)
signals = c(1,1.5)
pcs = c(30,50)
neighbours = c(10,25)
resolutions = c(1,1.5)

for(fragment in fragments){
  for(pct in pcts){
    for(signal in signals){
      for(pc in pcs){
        for(neighbour in neighbours){
          for(resolution in resolutions){
            print(paste(fragment,pct,signal,pc,neighbour,resolution, sep = "/"))
            scdata <- subset(
              x = smed_scatac,
              subset = peak_region_fragments < fragment & # tweak here? why removing cells with high number of peaks though?
                pct_reads_in_peaks > pct & # tweak here?
                nucleosome_signal < signal # tweak here?
            )
            scdata <- RunTFIDF(scdata)
            scdata <- FindTopFeatures(scdata, min.cutoff = 'q0')
            scdata <- RunSVD(scdata)
            scdata <- RunUMAP(object = scdata, reduction = 'lsi', dims = 2:pc)
            scdata <- FindNeighbors(object = scdata, reduction = 'lsi', dims = 2:pc, k.param = neighbour)
            scdata <- FindClusters(object = scdata, verbose = FALSE, algorithm = 3, resolution = resolution)
            umap <- DimPlot(object = scdata, label = TRUE) + ggtitle(paste(fragment,pct,signal,pc,neighbour,resolution, sep = " / "))
            pdf(
              file = 
                paste0(
                  "graphics/scatac_plots/iterations_cluster_scatac/",
                  paste("umap",fragment,pct,signal,pc,neighbour,resolution,sep = "_"),
                  ".pdf"
                  )
            )
            print(umap)
            dev.off()
            rm(scdata)
          }
        }
      }
    }
  }
}