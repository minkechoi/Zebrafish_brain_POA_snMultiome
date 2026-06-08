tgs_lit <- DF_tgs[!(DF_tgs$name %in% c("-",""," ")),]
tgs_lit$is <- 1

m <-
  tgs_lit[,c(8,9,10)] %>%
  pivot_wider(names_from = fate_tf, values_from = is) %>%
  column_to_rownames(var = "name") %>%
  as.matrix()

m[is.na(m)] <- 0

h_pre_abs <-
  Heatmap(
    name = "presence/absence",
    m,
    col = c("white","black"),
    cluster_columns = FALSE,
    top_annotation = 
      HeatmapAnnotation(
        fate = rep(names(broadcols_inf), each = 5),
        col = list(fate = broadcols_inf),
        show_legend = FALSE
      ),
    bottom_annotation = 
      HeatmapAnnotation(
        fate = rep(names(broadcols_inf), each = 5),
        col = list(fate = broadcols_inf),
        show_legend = FALSE
      ),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    column_names_side = "top"
  )

pdf("graphics/heatmap_presence_absence_targets_literature.pdf", height = 26, width = 10)
draw(h_pre_abs)
dev.off()














tgs_fig4 <- read.delim2("outputs/functional_annotation/targets_fig4_stripcharts.tsv")
DF_tgs_mainfig4 <- DF_tgs[DF_tgs$fate_tf %in% tgs_fig4$fate_tf,]
length(unique(DF_tgs_mainfig4$fate_tf))

ppp_mainfig4 <- list()
for(i in levels(DF_tgs_mainfig4$fate)){
  # i = "protonephridia"
  message(i)
  d <- DF_tgs_mainfig4[DF_tgs_mainfig4$fate == i, ]
  
  p_l <- list()
  for(j in unique(d$tf_sym)){
    # j = "foxA1"
    d_ = d[d$tf_sym == j,]
    
    # l_ <- length(unique(d_$name[!(d_$name %in% c(" ","-",""))]))
    # 
    # if(l_ > 30 ){
    #   set.seed(1234)
    #   margin_text <- rev(sort(
    #     sample(unique(d_$name[!(d_$name %in% c(" ","-",""))]),30)
    #   ))  
    # } else{
    #   margin_text <- rev(sort(unique(d_$name[!(d_$name %in% c(" ","-",""))])))
    # }
    
    margin_text <- d_$name[d_$name %in% tgs_fig4$tg_name[tgs_fig4$fate_tf == unique(d_$fate_tf)] ]
    
    margin_text <- margin_text[match(d_$name[order(d_$value,decreasing = TRUE)], margin_text)]
    margin_text <- margin_text[complete.cases(margin_text)]
    
    d_$in_text <- ifelse(d_$name %in% margin_text, TRUE, FALSE)
    
    set.seed(1234)
    p_l[[j]] <- ggplot(d_ %>% arrange(in_text), aes(x = tf_sym, y = value, fill = color, color = in_text, size = in_text)) +
      geom_jitter(width = 0.3, alpha = 0.7, pch = 21) +
      scale_size_manual(values = c(1.5,2.5))+
      scale_fill_identity() +
      scale_color_manual(values = c("white","black"))+
      labs(x = "TF Symbol", y = "Score (top 5% targets)") +
      theme_minimal() +
      theme(plot.margin = unit(c(0, 0, 0, 0), "in"))+
      guides(color ="none", size = "none")+
      grid.text(
        label = margin_text, x = .4, y = seq(0.20, 0.95, length.out = length(margin_text)), 
        gp = gpar(fontsize = 6,fontface = "italic", hjust = 1), draw = FALSE
      )
    
    
  }
  
  ppp_mainfig4[[i]] <- cowplot::plot_grid(plotlist = p_l, ncol = 2, nrow = 1)
}

pdf("graphics/ananse_top_stripchart_mainfig4.pdf",height = 3.5,width = 3.5)
for(i in names(ppp_mainfig4)){
  print(ppp_mainfig4[[i]])
}
dev.off()
















# TF CONNECTIVITY
tfs_kme_intramodular <- smed_id_module[smed_id_module$id %in% smed_tfs$id,]
tfs_kme_intramodular$kME <- 0

tf_eigen_ <- 
  WGCNA::signedKME(
    scale(t(smed_tfs_cw)), # all tfs, not only those with CV > 1.25 as in wgcna markdown
    MEs, outputColumnName = ""
  )

for(i in tfs_kme_intramodular$id){
  which_col <- which(colnames(tf_eigen_) == as.character(tfs_kme_intramodular$module[tfs_kme_intramodular$id == i]))
  tfs_kme_intramodular[tfs_kme_intramodular$id == i,3] <- tf_eigen_[i, which_col]
}


# TF CENTRALITY
smed_graph0_parse <- ParseNetwork(smed_graph0,list_attr = smed_attributes_list)

smed_id_module_ <- smed_id_module
colnames(smed_id_module_) <- c("id","member")

smed_grns_list <- divide_into_components( 
  x = smed_graph0_parse[[1]],
  CCs = smed_id_module_
)

smed_grns_list <-
  lapply(
    smed_grns_list,
    function(x){
      y = subgraph.edges(x,which(E(x)$weight > quantile(E(x)$weight,0.75)),delete.vertices = TRUE)
      return(y)
    }
  )

tfscentr_by_module <- centrality_by_network(smed_grns_list, normalized = TRUE)

top_central_tfs <- 
  unlist(
    sapply(
      tfscentr_by_module,
      function(x){
        if(length(x) < 3) {names(x)} else {names(rev(sort(x))[1:3])}
      }
    )
  )

smed_tfs_centrality_df <- 
  merge(
    data.frame(
      id = sub(".*\\.","",names(unlist(tfscentr_by_module))),
      module = sub("\\..*","",names(unlist(tfscentr_by_module))),
      centrality = unlist(tfscentr_by_module)
    ),
    smed_modules_table,
    by.x = 2,
    by.y = 4,
    all.x = TRUE
  ) %>% 
  mutate(top_central = ifelse(id %in% top_central_tfs,TRUE,FALSE)) %>%
  remove_rownames

smed_tfs_centrality_df <- smed_tfs_centrality_df[,c(2,3,1,9,10,11)]

colnames(smed_tfs_centrality_df) <- c("id","centrality","module","color","cellcolor","top_central")

smed_tfs_centrality_df = smed_tfs_centrality_df[smed_tfs_centrality_df$module %in% names( table(smed_tfs_centrality_df$module)[ table(smed_tfs_centrality_df$module)>2]),]


# CORRELATION OF BOTH
kme_centr_all <- merge(
  tfs_kme_intramodular,
  smed_tfs_centrality_df[,-3],
  by.x = 1,
  by.y = 1,
  all = TRUE
)

kme_centr <- kme_centr_all[complete.cases(kme_centr_all),]
kme_centr <- 
  kme_centr %>% 
  group_by(module) %>% 
  mutate(rel_centr = relativise(centrality)) %>%
  mutate(rel_kme = relativise(kME))
kme_centr <- kme_centr[complete.cases(kme_centr),]

# all together plot
pdf("graphics/scatter_cor_kme_all.pdf", wi = 5, he = 5)
plot(
  kme_centr$rel_kme,
  kme_centr$rel_centr,
  main = "agreement connectivity/centrality",
  bty = "n",
  xlab = "relative intra-modular connectivity",
  ylab = "relative centrality in module graph",
  pch = 21, bg = kme_centr$color
)
abline(a=0,b=1,lty = 2, col = "gray",lwd = 1.2)
dev.off()


kME_cent_PCC_all <- cor(kme_centr$rel_kme,kme_centr$rel_centr)
kME_cent_PCC_module <- 
  sapply(
    as.character(unique(kme_centr$module)),
    function(x){
      cor(
        kme_centr$kME[kme_centr$module == x],
        relativise(kme_centr$centrality[kme_centr$module == x])
      )
    }
  )

## module-wise plot
pdf("graphics/scatter_cor_kme_modulewise.pdf", he = 15, wi = 15)
par(mfrow = c(6,6))
for(i in sort(unique(kme_centr$module))){
  plot(
    kme_centr$rel_kme[kme_centr$module == i],
    kme_centr$rel_centr[kme_centr$module == i],
    pch = 21,
    bg = unique(kme_centr$color[kme_centr$module == i]),
    main = i,
    bty = "n",
    xlab = "relative intra-modular connectivity",
    ylab = "relative centrality in module graph"
  )
}
par(mfrow = c(1,1))
dev.off()

# All combined
d <- data.frame(cor = c(kME_cent_PCC_all,kME_cent_PCC_module))
d$module = rownames(d)
d$col = translate_ids(d$module, smed_modules_table[,c(4,7)])
d$col[d$module == ""] = "grey80"
d$how = "module-wise"
d$how[d$module == ""] = "whole network"
d$how <- factor(d$how, levels = c("whole network", "module-wise"))

pdf("graphics/cor_kme_centr_jitter.pdf", he = 5, wi = 5)
set.seed(1234)
d %>%
  ggplot(aes(x = how, y = cor, col = module, fill = module))+
  geom_jitter(cex = 2.2,width = .2, pch = 21)+
  scale_color_manual(values = darken(d$col,.5))+
  scale_fill_manual(values = d$col)+
  theme_classic()+
  ggtitle("Correlation between intramodular connectivity and centrality")+
  theme(legend.position = "none")
dev.off()


tfs_inf_stripchart <- read.delim2(
  file = "outputs/functional_annotation/tfs_inf_stripchart.tsv", header = T
)
list_inf_dfs <- lapply(lg_inf, function(x) x$influence)  # Extract data frames
# Create a vector of original sublist names corresponding to each row in data frames
fates <- unlist(lapply(seq_along(list_inf_dfs), function(i) rep(names(lg_inf)[i], nrow(list_inf_dfs[[i]]))))
inf_df_all <- do.call(rbind, list_inf_dfs)
inf_df_all$fate <- fates
inf_df_all$genesymbol <- 
  translate_ids(x = inf_df_all$factor, dict = tfs_all[,c(1,11)])

inf_df_top <- inf_df_top %>% group_by(fate) %>% slice_max(order_by = influence_score, n = 5)

list_tgs <- list()
for(i in unique(inf_df_top$fate)){
  
  fate <- unlist(strsplit(i, "_"))[1]
  
  tfs <- inf_df_top$factor[inf_df_top$fate == i]
  
  list_tgs[[fate]] <- 
    do.call(
      "rbind",
      lapply(
        tfs,
        function(x){
          y = E(lg[[fate]])[.from(x)]$prob # values
          names(y) = head_of(lg[[fate]], es = E(lg[[fate]])[.from(x)])$name # names
          y = y[y > quantile(y, .95)] # filter
          z = data.frame( # make a DF out of this
            fate = fate,
            tf = x,
            gene = names(y),
            value = y,
            value_rel = relativise(y)
          )
          rownames(z) = NULL
          return(z)
        }
      )
    )
}

DF_tgs <- do.call("rbind",list_tgs) # genius

DF_tgs$tf_sym <- translate_ids(x=DF_tgs$tf, dict = tfs_all[,c(1,11)])
DF_tgs$tf_sym[DF_tgs$tf_sym == "-"] <- DF_tgs$tf[DF_tgs$tf_sym == "-"]

DF_tgs$color <-
  broadcols_inf[
    match(
      DF_tgs$fate,
      names(broadcols_inf)
    )
  ]

DF_tgs$fate <- factor(DF_tgs$fate, levels = names(broadcols_inf))

DF_tgs <- DF_tgs[order(DF_tgs$fate, DF_tgs$tf),]

DF_tgs$fate_tf <- paste(as.character(DF_tgs$fate), DF_tgs$tf_sym)
DF_tgs$fate_tf <- factor(DF_tgs$fate_tf, levels = unique(DF_tgs$fate_tf))

DF_tgs$name <- translate_ids(x=DF_tgs$gene, dict = gene_names_lit,return.missing = FALSE)
DF_tgs$name[is.na(DF_tgs$name)] <- ""
head(DF_tgs)

ppp <- list()
for(i in levels(DF_tgs$fate)){
  message(i)
  d <- DF_tgs[DF_tgs$fate == i, ]
  
  p_l <- list()
  for(j in unique(d$tf_sym)){
    d_ = d[d$tf_sym == j,]
    
    l_ <- length(unique(d_$name[!(d_$name %in% c(" ","-",""))]))
    
    if(l_ > 30 ){
      set.seed(1234)
      margin_text <- rev(sort(
        sample(unique(d_$name[!(d_$name %in% c(" ","-",""))]),30)
      ))  
    } else{
      margin_text <- rev(sort(unique(d_$name[!(d_$name %in% c(" ","-",""))])))
    }
    
    margin_text <- margin_text[match(d_$name[order(d_$value,decreasing = TRUE)], margin_text)]
    margin_text <- margin_text[complete.cases(margin_text)]
    
    d_$in_text <- ifelse(d_$name %in% margin_text, TRUE, FALSE)
    
    d_ratio_ <- paste0(round(l_ / length(unique(DF_tgs$name)),2)*100,"%")
    
    set.seed(1234)
    p_l[[j]] <- ggplot(d_ %>% arrange(in_text), aes(x = tf_sym, y = value, fill = color, color = in_text, size = in_text)) +
      geom_jitter(width = 0.3, alpha = 0.7, pch = 21) +
      scale_size_manual(values = c(1.5,2.5))+
      scale_fill_identity() +
      scale_color_manual(values = c("white","black"))+
      labs(x = "TF Symbol", y = paste0("Score (top 5% targets) (",d_ratio_," in lit.)")) +
      theme_minimal() +
      theme(plot.margin = unit(c(0, 0, 0, 0), "in"))+
      guides(color ="none", size = "none")+
      grid.text(
        label = margin_text, x = .4, y = seq(0.20, 0.95, length.out = length(margin_text)), 
        gp = gpar(fontsize = 6,fontface = "italic", hjust = 1), draw = FALSE
      )
    
    
  }
  
  ppp[[i]] <- plot_grid(plotlist = p_l, ncol = 5, nrow = 1)
}

pdf("",height = 3.5,width = 12)
for(i in names(ppp)){
  print(ppp[[i]])
}
dev.off()



rosetta_ext_for_paper =
  rosetta_ext

rosetta_ext_for_paper$Symbol[rosetta_ext_for_paper$gene == "h1SMcG0016896"] = "foxf-1"
rosetta_ext_for_paper$name[rosetta_ext_for_paper$gene == "h1SMcG0016896"] = "foxf-1"
rosetta_ext_for_paper$name_lit[rosetta_ext_for_paper$gene == "h1SMcG0016896"] = "foxf-1"

rosetta_ext_for_paper =
  merge(
    rosetta_ext_for_paper,
    smed_id_module,agat_convert_sp_gff2gtf.pl
    by.x = 1,
    by.y = 1,
    all.x = TRUE
  )

write.table(
  rosetta_ext_for_paper,
  "outputs/20240605_supplementary_file_1_1.tsv",
  sep = "\t",
  dec = ".",
  row.names = TRUE,
  quote = FALSE
)

