
####FIG5D

DE_DF_tgs_mt_t3=list()

for(i in 1:6){
  dt=DE_DF_tgs[which(DE_DF_tgs$comp == unique(DE_DF_tgs$comp)[i]),]
  tf_d=unique(dt$tf)
  tg_d=unique(dt$tg)
  
  mt=matrix(0, nrow = length(tf_d),ncol = length(tg_d))
  rownames(mt)=tf_d
  colnames(mt)=tg_d
  
  for (k in 1:nrow(dt)) {
    mt[dt$tf[k],dt$tg[k]]=dt$value[k]
  }
  DE_DF_tgs_mt_t3[[unique(DE_DF_tgs$comp)[i]]]=mt

}


m_DE_DF_tgs_mt_t3=list()
p2=list("a"=c(1,2),"b"=c(1,2),
        "c"=c(3,4),"d"=c(3,4),
        "e"=c(5,6),"f"=c(5,6)
        )

for(i in 1:6){
  dt=DE_DF_tgs[which(DE_DF_tgs$comp == unique(DE_DF_tgs$comp)[i]),]
  comb_mt_tf = unique(c(rownames(DE_DF_tgs_mt_t3[[p2[[i]][1]]]),rownames(DE_DF_tgs_mt_t3[[p2[[i]][2]]]) ))

mt=matrix(0, nrow = length(comb_mt_tf),ncol = length(unique(dt$tg)))
rownames(mt)=comb_mt_tf
colnames(mt)=unique(dt$tg)

for (k in 1:nrow(dt)) {
  mt[dt$tf[k],dt$tg[k]]=dt$value[k]
}

m_DE_DF_tgs_mt_t3[[unique(DE_DF_tgs$comp)[i]]]=mt

}


library(openxlsx)

write.xlsx(
  m_DE_DF_tgs_mt_t3,
  file = "./outputs/danio_ANANSE_DEG_TF_mtrix.xlsx",
  #sheetName = comps,
  colNames = T, rowNames = TRUE, showNA = TRUE
)


##
ref_colors <- c("#0e6655", "#229954", "#138d75")

lighten_steps <- function(color, n = 2) {
  lighten_seq <- seq(0, 0.6, length.out = n)  
  lighten(color, amount = lighten_seq)
}

color_list <- lapply(ref_colors, lighten_steps)

all_colors <- unlist(color_list)

pdf(paste0("./figures/DE_DF_tgs_mt_hm_merged.pdf"), width =11, height = 9)
for (i in 0:2) {
set.seed(5678)

  a= Heatmap(
    name = unique(DE_DF_tgs$comp)[(2*i)+1],
    m_DE_DF_tgs_mt_t3[[(2*i)+1]],
    col = sequential_hcl(10, "Sunset"), #viridis(10),#colorRampPalette(c("#f1f5ff","#b1b8c4","#2c2b46"))(10), # sequential_hcl(10, "Sunset"),
    show_row_names = T,
    show_column_names = T, column_names_side = "bottom",
    column_km =2 ,
    column_names_gp = gpar(cex = .7),
    cluster_rows = T,
    clustering_method_rows = "average",
    #clustering_method_columns = "average", #"ward.D2",
    row_names_gp = gpar(cex = .7),
    #row_title_rot = 0,
    show_row_dend = T,
    km=2,
    row_title_gp = gpar(cex=.7)
  )
  
  b= Heatmap(
    name = unique(DE_DF_tgs$comp)[(2*i)+2],
    m_DE_DF_tgs_mt_t3[[(2*i)+2]],
    col = sequential_hcl(10, "Sunset"),#viridis(10),#colorRampPalette(c("#f1f5ff","#b1b8c4","#2c2b46"))(10), # sequential_hcl(10, "Sunset"),
    show_row_names = T,
    show_column_names = T, column_names_side = "bottom",
    column_km = 2,
    column_names_gp = gpar(cex = .7),
    #cluster_rows = T,
    #clustering_method_rows = "complete",
    #clustering_method_columns = "average", #"ward.D2",
    row_names_gp = gpar(cex = .7),
    #row_title_rot = 0,
    show_row_dend = T,
    km=2,
    row_title_gp = gpar(cex=.7)
  )
  
draw(a+b)
}
dev.off()





####fig5f

# Install if needed
# install.packages("igraph")


library(igraph)

angle <- -90 * pi / 180 # Convert 90 degrees to radians
# Rotation matrix for 2D (counter-clockwise)
rot_matrix <- matrix(c(cos(angle), sin(angle),
                       -sin(angle), cos(angle)), ncol = 2)

## ---- Node and edge definitions for cluster 35.0 ----

nodes_35.0 <- data.frame(
  id = c(
    "Stress","GC","cAMP","Ca2","MAPK","NR3C1",
    "CREB","CRTC_EP300","MEF2","KLF","PBX","POU2F2A","CHD",
    "NPAS4A_TF","FOS","EGR_TF","JUN","NR4A1_TF",
    "adgrb1b","fkbp5","npas4a_eff","plk2b","per2",
    "zgc122979","snrkb","egr1_eff","btg2","egr4_eff","si_dkey",
    # bPAC-specific effectors
    "nr1d2a","csrnp1b","diras1a","sat1a2","spsb4a","tp53bp2b","irs2b"
  ),
  label = c(
    "Acute stress","Glucocorticoids","cAMP/PKA","Ca2+ influx","MAPK/ERK","nr3c1 (GR)",
    "CREB/CREM","CRTC/EP300","MEF2","KLFs","PBX","pou2f2a","CHD1/2",
    "npas4a (TF)","FOS/FOSL2","EGR1/EGR4 (TF)","JUN","nr4a1 (TF)",
    "adgrb1b","fkbp5","npas4a (eff)","plk2b","per2",
    "zgc:122979","snrkb","egr1 (eff)","btg2","egr4 (eff)","si:dkey-27j5.5",
    "nr1d2a","csrnp1b","diras1a","sat1a.2","spsb4a","tp53bp2b","irs2b"
  ),
  group = c(
    rep("shared", 29),          # first 29 nodes shared
    rep("bPAC", 7)              # last 7 nodes bPAC-specific
  ),
  stringsAsFactors = FALSE
)

edges_35.0 <- data.frame(
  from = c(
    "Stress","Stress","Stress",
    "Ca2","GC",
    "cAMP","cAMP","CRTC_EP300","Ca2",
    "MAPK","MAPK","MAPK",
    "NR3C1","NR3C1","NR3C1",
    "CREB","CREB","CREB",
    "MEF2","KLF","PBX","POU2F2A","CHD",
    "NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF",
    "FOS","FOS","FOS",
    "EGR_TF","EGR_TF","EGR_TF",
    "NR3C1","NR3C1","NR3C1","NR3C1","NR3C1",
    "NPAS4A_TF","NPAS4A_TF","NPAS4A_TF"
  ),
  to = c(
    "GC","cAMP","Ca2",
    "MAPK","NR3C1",
    "CREB","CRTC_EP300","CREB","CREB",
    "FOS","EGR_TF","JUN",
    "NPAS4A_TF","FOS","EGR_TF",
    "NPAS4A_TF","FOS","EGR_TF",
    "NPAS4A_TF","EGR_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF",
    "btg2","npas4a_eff","plk2b","per2","zgc122979","snrkb","egr4_eff","si_dkey",
    "btg2","plk2b","adgrb1b",
    "btg2","plk2b","per2",
    "fkbp5","nr1d2a","csrnp1b","irs2b","tp53bp2b",
    "sat1a2","diras1a","spsb4a"
  ),
  stringsAsFactors = FALSE
)

## ---- Build graph and set styles ----

g35.0 <- graph_from_data_frame(
  d = edges_35.0,
  vertices = nodes_35.0,
  directed = TRUE
)

# Colors: pink for bPAC-specific nodes, grey for shared
V(g35.0)$color <- ifelse(V(g35.0)$group == "bPAC", "pink", "grey80")
V(g35.0)$label <- V(g35.0)$label
V(g35.0)$size  <- ifelse(V(g35.0)$group == "bPAC", 20, 15)
E(g35.0)$arrow.size <- 0.4

# Layout: layered tree-ish layout from upstream to downstream
# You can manually adjust or use layout_as_tree with a chosen root.
set.seed(123)
lay_35.0 <- layout_with_sugiyama(g35.0)$layout
#lay_35.0 <- lay_35.0 %*% rot_matrix # Apply rotation

## ---- Export to PDF ----

pdf("./figures/cluster_35.0_pathway.pdf", width = 8, height = 6)
plot(
  g35.0,
  layout = lay_35.0,
  vertex.label.cex = 0.7,
  vertex.label.color = "black",
  main = "Cluster 35.0 acute-stress pathway"
)
legend(
  "topleft",
  legend = c("Shared node", "bPAC-specific node"),
  col = c("grey80", "pink"),
  pch = 19,
  pt.cex = 1,
  bty = "n"
)
dev.off()



## ---- Node and edge definitions for cluster 35.1 ----

nodes_35.1 <- data.frame(
  id = c(
    "Stress","GC","cAMP","Ca2","MAPK","NR3C1",
    "CREB","CRTC_EP300","MEF2","KLF","PBX","POU2F2A","CHD",
    "NPAS4A_TF","FOS","EGR_TF","JUN","NR4A1_TF",
    # Shared effectors (present in control & bPAC)
    "itm2cb","usp2a","egr4_eff","btg2","plk3","crema_eff",
    "atp1b1b","fosaa_eff","fosl2_eff","hip1",
    # bPAC-specific effectors
    "fkbp5","npas4a_eff","per2","nr4a1_eff","nr1i2",
    "sgsm3","plk2b","sat1a2","ptp4a2b","smad7","spry4","angpt1",
    "drd4","nts","pcsk1","irs2b","scg3"
  ),
  label = c(
    "Acute stress","Glucocorticoids","cAMP/PKA","Ca2+ influx","MAPK/ERK","nr3c1 (GR)",
    "CREB/CREM","CRTC/EP300","MEF2","KLFs","PBX","pou2f2a","CHD1/2",
    "npas4a (TF)","FOS/FOSL2","EGR1/EGR4 (TF)","JUN","nr4a1 (TF)",
    "itm2cb","usp2a","egr4 (eff)","btg2","plk3","crema (eff)",
    "atp1b1b","fosaa (eff)","fosl2 (eff)","hip1",
    "fkbp5","npas4a (eff)","per2","nr4a1 (eff)","nr1i2",
    "sgsm3","plk2b","sat1a.2","ptp4a2b","smad7","spry4","angpt1",
    "drd4-rs","nts","pcsk1","irs2b","scg3"
  ),
  group = c(
    rep("shared", 29),
    rep("bPAC", 16)
  ),
  stringsAsFactors = FALSE
)

edges_35.1 <- data.frame(
  from = c(
    "Stress","Stress","Stress",
    "Ca2","GC",
    "cAMP","cAMP","CRTC_EP300","Ca2",
    "MAPK","MAPK","MAPK",
    "NR3C1","NR3C1","NR3C1",
    "CREB","CREB","CREB",
    "MEF2","KLF","PBX","POU2F2A","CHD",
    # shared effector edges
    "NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF",
    "NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF",
    "FOS","FOS","FOS","FOS",
    "EGR_TF","EGR_TF",
    # bPAC GR/CREB outputs
    "NR3C1","NR3C1","NR3C1",
    "CREB","CREB","CREB",
    # bPAC effector edges from NPAS4A/FOS/EGR
    "NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF",
    "FOS","FOS","FOS","FOS",
    "NPAS4A_TF"
  ),
  to = c(
    "GC","cAMP","Ca2",
    "MAPK","NR3C1",
    "CREB","CRTC_EP300","CREB","CREB",
    "FOS","EGR_TF","JUN",
    "NPAS4A_TF","FOS","EGR_TF",
    "NPAS4A_TF","FOS","EGR_TF",
    "NPAS4A_TF","EGR_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF",
    "btg2","plk3","crema_eff","itm2cb","atp1b1b",
    "fosaa_eff","fosl2_eff","hip1","egr4_eff",
    "btg2","plk3","fosaa_eff","fosl2_eff",
    "egr4_eff","ptp4a2b",
    "fkbp5","per2","nr1i2",
    "fkbp5","per2","nr1i2",
    "sgsm3","plk2b","sat1a2","ptp4a2b","spry4","angpt1","irs2b","npas4a_eff",
    "drd4","nts","pcsk1","ptp4a2b",
    "scg3"
  ),
  stringsAsFactors = FALSE
)

## ---- Build graph and style ----

g35.1 <- graph_from_data_frame(
  d = edges_35.1,
  vertices = nodes_35.1,
  directed = TRUE
)

V(g35.1)$color <- ifelse(V(g35.1)$group == "bPAC", "pink", "grey80")
V(g35.1)$label <- V(g35.1)$label
V(g35.1)$size  <- ifelse(V(g35.1)$group == "bPAC", 18, 15)
E(g35.1)$arrow.size <- 0.4

set.seed(123)

lay_35.1 <- layout_with_sugiyama(g35.1)$layout
#lay_35.1 <- lay_35.1 %*% rot_matrix # Apply rotation


## ---- Export to PDF ----



pdf("./figures/cluster_35.1_pathway.pdf", width = 8, height = 6)
plot(
  g35.1,
  layout = lay_35.1,
  vertex.label.cex = 0.7,
  vertex.label.color = "black",
  main = "Cluster 35.1 acute-stress pathway"
)
legend(
  "topleft",
  legend = c("Shared node", "bPAC-specific node"),
  col = c("grey80", "pink"),
  pch = 19,
  pt.cex = 1.5,
  bty = "n"
)
dev.off()



## ---- Node and edge definitions for cluster 45 ----

nodes_45 <- data.frame(
  id = c(
    "Stress","GC","cAMP","Ca2","MAPK","NR3C1",
    "CREB","CRTC_EP300","MEF2","KLF","PBX","POU2F2A","CHD",
    "NPAS4A_TF","FOS","EGR_TF","JUN","NR4A1_TF",
    # shared effectors
    "per2","btg2","adgrb1b","npas4a_eff","fosl2_eff",
    # bPAC-specific effectors
    "fkbp5","galr2b","ptp4a2b",
    "slc20a1b","slc23a2","slc25a25b","slc38a3a","slc7a14a",
    "scg3","arg2","smox","sat1a2","sgsm3",
    "sik1","sik2b","skilb","smad7","trib2","ptges3a",
    "bmp1a","kcnk10a","nfil3","usp2a","zfand5b"
  ),
  label = c(
    "Acute stress","Glucocorticoids","cAMP/PKA","Ca2+ influx","MAPK/ERK","nr3c1 (GR)",
    "CREB/CREM","CRTC/EP300","MEF2","KLFs","PBX","pou2f2a","CHD1/2",
    "npas4a (TF)","FOS/FOSL2","EGR1/EGR4 (TF)","JUN","nr4a1 (TF)",
    "per2","btg2","adgrb1b","npas4a (eff)","fosl2 (eff)",
    "fkbp5","galr2b","ptp4a2b",
    "slc20a1b","slc23a2","slc25a25b","slc38a3a","slc7a14a",
    "scg3","arg2","smox","sat1a.2","sgsm3",
    "sik1","sik2b","skilb","smad7","trib2","ptges3a",
    "bmp1a","kcnk10a","nfil3-6","usp2a","zfand5b"
  ),
  group = c(
    rep("shared", 23),
    rep("bPAC", 24)
  ),
  stringsAsFactors = FALSE
)

edges_45 <- data.frame(
  from = c(
    "Stress","Stress","Stress",
    "Ca2","GC",
    "cAMP","cAMP","CRTC_EP300","Ca2",
    "MAPK","MAPK","MAPK",
    "NR3C1","NR3C1","NR3C1",
    "CREB","CREB","CREB",
    "MEF2","KLF","PBX","POU2F2A","CHD",
    # shared effector edges
    "NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF",
    "FOS","FOS","EGR_TF",
    "NR3C1","NR3C1","NR3C1",
    # bPAC-specific effector edges
    "FOS","NPAS4A_TF",
    "FOS","FOS","FOS","FOS","FOS",
    "NPAS4A_TF","NPAS4A_TF","EGR_TF",
    "NPAS4A_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF",
    "FOS","FOS",
    "NPAS4A_TF","NPAS4A_TF",
    "FOS","FOS","NPAS4A_TF","FOS","NPAS4A_TF","EGR_TF",
  ),
  to = c(
    "GC","cAMP","Ca2",
    "MAPK","NR3C1",
    "CREB","CRTC_EP300","CREB",
    "FOS","EGR_TF","JUN",
    "NPAS4A_TF","FOS","EGR_TF",
    "NPAS4A_TF","FOS","EGR_TF",
    "NPAS4A_TF","EGR_TF","NPAS4A_TF","NPAS4A_TF","NPAS4A_TF",
    "btg2","npas4a_eff","per2","fosl2_eff",
    "btg2","adgrb1b","per2",
    "per2","fkbp5","smad7",
    "galr2b","galr2b",
    "slc20a1b","slc23a2","slc25a25b","slc38a3a","slc7a14a",
    "scg3","sgsm3","sat1a2",
    "sik1","sik2b","skilb","nfil3",
    "arg2","smox",
    "trib2","ptges3a",
    "bmp1a","kcnk10a","usp2a","zfand5b","ptp4a2b","ptp4a2b","ptp4a2b"
  ),
  stringsAsFactors = FALSE
)

## ---- Build graph and style ----

g45 <- graph_from_data_frame(
  d = edges_45,
  vertices = nodes_45,
  directed = TRUE
)

V(g45)$color <- ifelse(V(g45)$group == "bPAC", "pink", "grey80")
V(g45)$label <- V(g45)$label
V(g45)$size  <- ifelse(V(g45)$group == "bPAC", 18, 15)
E(g45)$arrow.size <- 0.4

set.seed(123)
lay_45 <- layout_with_sugiyama(g45)$layout
#lay_45 <- lay_45 %*% rot_matrix # Apply rotation

## ---- Export to PDF ----

pdf("./figures/cluster_45_pathway.pdf", width = 8, height = 6)
plot(
  g45,
  layout = lay_45,
  vertex.label.cex = 0.7,
  vertex.label.color = "black",
  main = "Cluster 45 acute-stress pathway"
)
legend(
  "topleft",
  legend = c("Shared node", "bPAC-specific node"),
  col = c("grey80", "pink"),
  pch = 19,
  pt.cex = 1.5,
  bty = "n"
)
dev.off()
