library(presto)
library(GenomicRanges)
library(viridis)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)

set.seed(1234)
options(future.globals.maxSize = 120000 * 1024^2)
r.variable=4000
vs="06_2025_v4"
umap="wnn.umap"

#set working dir
if(getwd()!=paste0("D:/projects/scMultiome_oxt/",vs)){
  setwd(paste0("./",vs))
}

#functions adopted and modified from Alberto Perez-Posada @apposada
source("./ext_code/r_code/functions/sourcefolder.R")
sourceFolder("./ext_code/r_code/functions/")


ss1= readRDS("./data/rds/Step8_motif_var4000.rds")
#bPAC

  #find top DARs
ids=list("bPAC_post45_vs_cont_pre"=c("45_sst1.1:bPAC_post", "45_sst1.1:cont_pre"),
         "bPAC_post45_vs_cont_post"=c("45_sst1.1:bPAC_post", "45_sst1.1:cont_post"),
         "bPAC_post45_vs_bPAC_pre"=c("45_sst1.1:bPAC_post", "45_sst1.1:bPAC_pre"),
         "bPAC_post351_vs_cont_pre"=c("35.1_avp:bPAC_post", "35.1_avp:cont_pre"),
         "bPAC_post351_vs_cont_post"=c("35.1_avp:bPAC_post", "35.1_avp:cont_post"),
         "bPAC_post351_vs_bPAC_pre"=c("35.1_avp:bPAC_post", "35.1_avp:bPAC_pre")
         )

dapeak_list=list()
for (i in 1:length(ids)) {
  da_peaks <- FindMarkers(
    object = ss1,
    assay = "ATAC_macs3",
    group.by = "merged_sub.anno_type_ori",
    ident.1 = ids[[i]][1],
    ident.2 = ids[[i]][2],
    only.pos = TRUE,
    test.use = 'LR',
    min.pct = 0.05,
    latent.vars = 'nCount_ATAC_macs3'
  )
  
  dar=da_peaks
  top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.05 & da_peaks$pct.1 > 0.02, ])  
  top.da.peak.loc <- da_peaks[da_peaks$p_val < 0.05 & da_peaks$pct.1 >  0.02, ]
  top.peak.links=as.data.frame(ss1@assays$ATAC_macs3@links[which(ss1@assays$ATAC_macs3@links$peak %in% top.da.peak )])
  top.peak_genes=unique(top.peak.links$gene)
  
  dapeak_list[[i]]=list("dar"=dar,"peak"=top.da.peak,"peak_all"=top.da.peak.loc,"peak_anno"=top.peak.links,"genes"=top.peak_genes)
}
  save(
  dapeak_list,
  file = "./data/rda/top_peak_45.rda"
)

 ' result=list(
  "dar"=da_peaks,  
  "top.da.peak" <- top.da.peak,  
  "top.da.peak.loc" <- da_peaks[da_peaks$p_val < peak_p_val & da_peaks$pct.1 > pct1, ],
  "top.peak.links"=ss1@assays$assay_type@links[which(ss1x@assays$assay_type@links$peak %in% top.da.peak )],
  "top.peak_genes"=unique(top.peak.links$gene)
  )'



# get top differentially accessible peaks

deg45_bPAC=read.csv("./outputs/DEGs/DEGs_postLD_preLD_45.sst1.1.bPAC.post_45.sst1.1.bPAC.pre.csv",row.names = 1)
deg45_bPAC=dplyr::filter(deg45_bPAC, padj <0.05)

deg45_post=read.csv("./outputs/DEGs/DEGs_bPAC_cont_45.sst1.1.bPAC.post_45.sst1.1.cont.post.csv",row.names = 1)
deg45_post=dplyr::filter(deg45_post, pvalue <0.005)


intersect(sst.top.peak, rownames(deg45_bPAC))

#bPAC
enriched.motifs <- FindMotifs(
  object = ss1,
  features = top.da.peak,
  background = 40000,
  assay = "ATAC_macs3",
  verbose = TRUE,
  p.adjust.method = "BH"
  
)

enriched.motifs=dplyr::filter(enriched.motifs,p.adjust <0.05)

a=MotifPlot(
  object = ss1,
  motifs = head(rownames(enriched.motifs),20)
)

#control
enriched.motifs_cont <- FindMotifs(
  object = ss1,
  features = top.da.peak_cont,
  background = 40000,
  assay = "ATAC_macs3",
  verbose = TRUE,
  p.adjust.method = "BH"
  
)

enriched.motifs_cont=dplyr::filter(enriched.motifs_cont,p.adjust <0.05)

aa=MotifPlot(
  object = ss1,
  motifs = head(rownames(enriched.motifs_cont),20)
)

a/aa
library(ggrepel)

DMR_enriched_TF=left_join(enriched.motifs_cont,enriched.motifs, by="motif.name")
DMR_enriched_TF[is.na(DMR_enriched_TF)]=0
DMR_enriched_TF_tb= DMR_enriched_TF[,c(8,6,15,9,17)]
DMR_enriched_TF_tb[,4]=-log10(DMR_enriched_TF_tb[,4])
DMR_enriched_TF_tb[,5]=-log10(DMR_enriched_TF_tb[,5])
DMR_enriched_TF_tb$p.adjust.y[is.infinite(DMR_enriched_TF_tb$p.adjust.y)]=0
DMR_enriched_TF_tb["delta_FC"]=DMR_enriched_TF_tb$fold.enrichment.y-DMR_enriched_TF_tb$fold.enrichment.x

p <- ggplot(DMR_enriched_TF_tb, aes(p.adjust.x, p.adjust.y,label=motif.name,size=abs(delta_FC), color=fold.enrichment.y))

pdf(paste0("./figures/ATAC/top_peaks/45_sst_enriched_motif_comp.pdf"),
    width = 8,height =8)
p + geom_point()+geom_abline(slope = 1, intercept = 0,linetype=3)+geom_text_repel()+ylim(c(0,7))+xlim(c(0,7))+theme_classic()
dev.off()

DMR_enriched_TF_tb[,"ratio"]= DMR_enriched_TF_tb$p.adjust.y/DMR_enriched_TF_tb$p.adjust.x
DMR_enriched_TF_bPAc=DMR_enriched_TF_tb$motif.name[c(which(DMR_enriched_TF_tb$ratio >= 1))]


#expressed
all_marks=read.csv("./outputs/cell_type/all.sub.markers_var4000_merged_sub.anno.csv",row.names = 1)
c45_mark=dplyr::filter(all_marks, cluster == "45_sst1.1")

en.motif= unique(c(sapply(str_split(string = tolower(DMR_enriched_TF_bPAc),
                                    pattern = "\\(|\\:"), `[`, 1),na.omit(sapply(str_split(string = tolower(DMR_enriched_TF_bPAc), pattern = "\\(|\\:"), `[`, 3)) ))

overlapped= list()

for (i in en.motif) {
  overlapped[[i]]=c45_mark$gene[grepl(pattern = i, x = c45_mark$gene)]
}

exp_enriched_motifs=c("NR3C2","Ar","CREM","Foxo1","FOSL2::JUN(var.2)","FOSL2::JUND(var.2)","FOSL2::JUNB(var.2)","KLF9","SP4")

pdf(paste0("./figures/ATAC/top_peaks/45_sst_exp_enriched_motif.pdf"),
    width = 12,height = 6)
MotifPlot(
  object = ss1,
  motifs = head(rownames(enriched.motifs[which(enriched.motifs$motif.name %in% exp_enriched_motifs),]),10)
)
dev.off()

FeatureDimPlot(
  srt = ss1, 
  features = c("nr3c2",  "ar", "foxo1"),
  assay = "SCT",
  seed = 0, compare_features = F, label =F, label_repel = T,label_insitu = TRUE, 
  add_density = F,palette = "GdRd",bg_cutoff = 0.5, 
  reduction = "wnn.umap", 
  theme_use = "theme_blank",
  theme_args=theme(strip.text = element_text(size = 20,face = c("bold.italic")))
)

####IEGS
library(dplyr)
IEGs = c("fosaa","fosab","fosb","fosl2","itm2cb","egr1","egr2a","egr2b","egr3","ier2a")
load(  file = "./data/rda/top_peak_45.rda")

dar_IEGs= c()
top.da.peak=c()
for (i in 1:6) {
  dar_IEGs=unique(c(dar_IEGs,intersect(dapeak_list[[i]]$genes,IEGs)))
  top.da.peak=unique(c(top.da.peak,row.names(dapeak_list[[i]]$dar)))
}

peak_ann=dapeak_list[[1]]$peak_anno
IEG_top.peak=as.data.frame(ss1@assays$ATAC_macs3@links) %>% dplyr::filter(peak %in% top.da.peak) %>% dplyr::filter(gene %in% IEGs)


enriched.motifs_IEG <- FindMotifs(
  object = ss1,
  features = IEG_top.peak$peak,
  background = 40000,
  assay = "ATAC_macs3"
)
write.csv(enriched.motifs_IEG,"./outputs/ATAC_res/IEG_motif.csv")

enriched.motifs_IEG_ft=dplyr::filter(enriched.motifs_IEG,p.adjust <0.05)

pdf(paste0("./figures/ATAC/top_peaks/45.351_sst_3IEG_enriched_motif.pdf"),
    width = 12,height = 6)
MotifPlot(
  object = ss1,
  motifs = head(rownames(enriched.motifs_IEG_ft),20)
)
dev.off()

###exp markers
all_marks=read.csv("./outputs/cell_type/all.sub.markers_var4000_merged_sub.anno.csv",row.names = 1)
c45_mark=dplyr::filter(all_marks, cluster == "45_sst1.1")
c351_mark=dplyr::filter(all_marks, cluster == "35.1_avp")
en.motif= unique(c(sapply(str_split(string = tolower(enriched.motifs_IEG_ft$motif.name),
                                    pattern = "\\(|\\:"), `[`, 1),na.omit(sapply(str_split(string = tolower(enriched.motifs_IEG_ft$motif.name), pattern = "\\(|\\:"), `[`, 3)) ))

overlapped= c()

for (i in en.motif) {
  overlapped=c(overlapped,
               c45_mark$gene[grepl(pattern = i, x = c45_mark$gene)],
               c351_mark$gene[grepl(pattern = i, x = c351_mark$gene)])
}
