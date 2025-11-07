library(CoVarNet)

# load data
scr<-readRDS("./data/seurat.rds")

# processing
sce <- subset(scr, subset = Tissue != "PBMC" & !cellType2 %in% c("Unknown", "Platelets"))
meta<-sce@meta.data
rt_sp <- names(table(meta$orig.ident))[table(meta$orig.ident) >= 100]
meta <- meta[meta$orig.ident %in% rt_sp, ]
meta$Origin2<-factor(meta$Origin2, levels = c("LN_P","LN_N", "Normal", "Adjacent", "Tumor"))
meta_filter <- meta[,c('orig.ident','cellType','cellType3','Patient','Origin2')]
colnames(meta_filter) <- c('sampleID','majorCluster','subCluster','cohort','group')

# choose k value
mat_fq_raw<-freq_calculate(meta_filter)
mat_fq_norm<-freq_normalize(mat_fq_raw,normalize="minmax") 
res <- nmf(mat_fq_norm, rank = 2:20, method = "nsNMF", seed = rep(123456, 6), .options = "vp")
plot(res)

# run nmf
K=12
NMF_K12<-nmf(mat_fq_norm, K, method = "nsNMF", seed = rep(77, 6), nrun = 30)
colnames(basis(NMF_K12))=c("CM01","CM02","CM03","CM04","CM05","CM06","CM07","CM08","CM09","CM10","CM11","CM12")
rownames(coef(NMF_K12))=c("CM01","CM02","CM03","CM04","CM05","CM06","CM07","CM08","CM09","CM10","CM11","CM12")

# figure 2e
pdf("nmf_character_plot_tissue_distribution.pdf", width = 10, height = 8)
gr.distribution(NMF_K12,meta=meta_filter,group="group") 
dev.off()

# Extended Data Figure8a
pdf("nmf_character_plot_all.pdf", width = 10, height = 12)
gr.weight_all(NMF_K12)
dev.off()

# Extended Data Figure8b
pdf("nmf_character_plot_top15.pdf", width = 12, height = 8)
gr.weight_top(NMF_K12,num=15)
dev.off()








