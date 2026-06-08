## "code/r_code/r_scripts/wgcna_module_exploration.R"

dir <- "~/projects/smed_cisreg/"
setwd(dir)

library(topGO)
source("code/r_code/r_functions/sourcefolder.R")
sourceFolder("code/r_code/functions/")

load("outputs/rda/wgcna_module_exploration.rda") # contains datExpr, smed_id_module_wgcna, smed_cisreg_ctypes, and whatever is needed

# Prep

smed_wg_module <- 
  merge(
    t(scale(datExpr)),
    smed_id_module_wgcna,
    by.x = 0, by.y = 1,
    all.X = TRUE
  )
rownames(smed_wg_module) <- smed_wg_module[,1]
smed_wg_module[,1] <- NULL

smed_wg_avgexp <-
  aggregate(
    smed_wg_module[,c(1:ncol(smed_wg_module)-1)],
    by = list(smed_wg_module$module),
    FUN = mean
  )
rownames(smed_wg_avgexp) <- smed_wg_avgexp[,1]
smed_wg_avgexp <- smed_wg_avgexp[,-1]

# Barplots of all the gene modules
pdf(
  file = "graphics/wgcna_exploration_modules/smed_cisreg_wgcna_module_exploration_barplots.pdf",
  width = 8, height = 5
  )
par(
  mar=c(12,4,4,2)+0.1,
  xpd = TRUE
  )

for (i in 1:nrow(smed_wg_avgexp)) {
  module_i <- rownames(smed_wg_avgexp)[i]
  dat <- unlist(c((smed_wg_avgexp[i, 1:ncol(smed_wg_avgexp)])))
  
  barplot(
    height = dat, col = smed_ctypes$col,#[match(colnames(smed_wg_avgexp),smed_ctypes$ctype)],
    las=2, cex.names=0.7, ylab = "z-score",
    main = paste0("Smed module ",module_i),
    )
}
dev.off()

# GO term analysis
# Gene universe
gene_universe <- rownames(smed_cpm)

# Gene-GO mappings
smed_id_GO <- readMappings("./outputs/gene_annotation/smed_GOs.tsv")

# List of genes of interest
smed_wg_list <- split(smed_id_module_wgcna$id,smed_id_module_wgcna$module)

# GO term analysis wrapper
smed_wg_GO_all <- 
  getGOs(
    smed_wg_list,
    gene_universe= rownames(smed_cpm),
    gene2GO = smed_id_GO
    )

pdf("graphics/wgcna_exploration_modules/smed_cisreg_wgcna_module_exploration_GOs.pdf")
for (i in smed_wg_GO_all$GOplot) {print(i)}
dev.off()