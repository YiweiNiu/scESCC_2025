
# Tidyverse
library(tidyverse)

# ggplot
library(cowplot)
library(ggpubr)
library(ggcorrplot)
library(ggrepel)
theme_set(theme_cowplot())

# heatmap
library(pheatmap)

# color
library(ggsci)
library(RColorBrewer)

# single cell
library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)

# doublets
library(scds)

# soupx
library(SoupX)

source(here::here("code/preprocess.R"))
source(here::here("code/plot.R"))

#这个脚本用于清洗数据
#* 过滤低质量细胞
#* 过滤双细胞
#* 使用 SoupX 矫正 count
#* 使用矫正后的 count 构建 Seurat 对象，修改细胞名字，并且保存到磁盘

cellranger_dir = '/niuyw-usb-disk/Projects/scESCA/cellranger/cellranger_count_outs'

# calculate confounders
.confounder_cal <- function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.heat"]] <- PercentageFeatureSet(seurat_obj, features = filter_features(seurat_obj, heat_genes.hs))
  seurat_obj[["percent.dissociation"]] <- PercentageFeatureSet(seurat_obj, features = filter_features(seurat_obj, dissocitaion_genes.hs))
  RPS.genes <- grep(pattern = "^RPS", x = rownames(seurat_obj), value = TRUE)
  RPL.genes <- grep(pattern = "^RPL", x = rownames(seurat_obj), value = TRUE)
  ribo.genes = c(RPS.genes, RPL.genes)
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, features = ribo.genes)
  return (seurat_obj)
}

# build seurat obj
.build_seurat <- function(filtered_count_dir, sam) {
  # seurat
  seurat.data = Read10X(data.dir = filtered_count_dir)
  seurat.raw = CreateSeuratObject(counts = seurat.data)
  seurat.raw = .confounder_cal(seurat_obj = seurat.raw)
  seurat.raw$Source = str_replace(sam, '_', '.')
  return(seurat.raw)
}

# filter cells
.filter_cells <- function(seurat_obj){
  # filter low-quality cells
  reasons = lowQuality_filter(seurat_obj, gene_fixed_low = 500, gene_fixed_high = 3000)
  seurat = seurat_obj[,!reasons$discard]

  # filter doublets
  seurat = .doublet_scds(seurat)
  high_scds = seurat$high_scds
  print(table(high_scds))
  seurat = subset(seurat, subset = high_scds == FALSE)
  return(seurat)
}

# soupx
.soupx <- function(raw_count_dir, seurat_obj) {
  seurat = fast_cluster(seurat_obj, res = 3)
  tod = Read10X(data.dir = raw_count_dir)
  toc = GetAssayData(seurat, assay = 'RNA', slot = 'counts')
  sc = SoupChannel(tod, toc, keepDroplets = TRUE)

  DR <- FetchData(seurat,
                  vars = c("seurat_clusters", "UMAP_1", "UMAP_2")) %>%
    dplyr::select(RD1 = UMAP_1, RD2 = UMAP_2, Cluster = seurat_clusters)

  # profiling the soup
  sc = estimateSoup(sc)
  # add dimenstion-reduction
  sc = setDR(sc, DR)
  # add cluster info
  sc = setClusters(sc, DR$Cluster)
  return(list(sc = sc, DR = DR))
}

# use mp markers to estimate soup
sams = c('S1125_P1')
for (sam in sams) {
  data_dir = file.path(cellranger_dir, sam)
  raw_count_dir = file.path(data_dir, 'raw_feature_bc_matrix')
  filtered_count_dir = file.path(data_dir, 'filtered_feature_bc_matrix')
  # prepare for soupx
  seurat.raw = .build_seurat(filtered_count_dir, sam)
  seurat = .filter_cells(seurat.raw)
  seurat

  tmp = .soupx(raw_count_dir, seurat)
  sc = tmp$sc
  DR = tmp$DR
  #plotMarkerDistribution(sc)

  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = soupx_genes$mp))
  #plotMarkerMap(sc, geneSet = soupx_genes$mp, DR = DR, useToEst = useToEst)

  sc = calculateContaminationFraction(sc, list(IG = soupx_genes$mp), useToEst = useToEst)
  mean(sc$metaData$rho)
  outs = adjustCounts(sc)
  seurat <- CreateSeuratObject(counts = outs)
  # save
  save_path = file.path('/niuyw-usb-disk/Projects/scESCA/cellranger/soupx', sam)
  dir.create(save_path, showWarnings = FALSE)
  saveRDS(seurat, file = file.path(save_path, 'seurat.rds'))
}

# use fibroblasts markers to estimate soup
sams = c('S0819_N4')
for (sam in sams) {
  data_dir = file.path(cellranger_dir, sam)
  raw_count_dir = file.path(data_dir, 'raw_feature_bc_matrix')
  filtered_count_dir = file.path(data_dir, 'filtered_feature_bc_matrix')
  # prepare for soupx
  seurat.raw = .build_seurat(filtered_count_dir, sam)
  seurat = .filter_cells(seurat.raw)
  seurat

  tmp = .soupx(raw_count_dir, seurat)
  sc = tmp$sc
  DR = tmp$DR
  #plotMarkerDistribution(sc)

  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = soupx_genes$fibroblasts))
  #plotMarkerMap(sc, geneSet = soupx_genes$fibroblasts, DR = DR, useToEst = useToEst)

  sc = calculateContaminationFraction(sc, list(IG = soupx_genes$fibroblasts), useToEst = useToEst)
  mean(sc$metaData$rho)
  outs = adjustCounts(sc)
  seurat <- CreateSeuratObject(counts = outs)
  # save
  save_path = file.path('/niuyw-usb-disk/Projects/scESCA/cellranger/soupx', sam)
  dir.create(save_path, showWarnings = FALSE)
  saveRDS(seurat, file = file.path(save_path, 'seurat.rds'))
}

# use platelets markers to estimate soup
sams = c('S0619_P2', 'S0819_P2')
for (sam in sams) {
  data_dir = file.path(cellranger_dir, sam)
  raw_count_dir = file.path(data_dir, 'raw_feature_bc_matrix')
  filtered_count_dir = file.path(data_dir, 'filtered_feature_bc_matrix')
  # prepare for soupx
  seurat.raw = .build_seurat(filtered_count_dir, sam)
  seurat = .filter_cells(seurat.raw)
  seurat

  tmp = .soupx(raw_count_dir, seurat)
  sc = tmp$sc
  DR = tmp$DR
  #plotMarkerDistribution(sc)

  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = soupx_genes$platelet))
  #plotMarkerMap(sc, geneSet = soupx_genes$platelet, DR = DR, useToEst = useToEst)

  sc = calculateContaminationFraction(sc, list(IG = soupx_genes$platelet), useToEst = useToEst)
  mean(sc$metaData$rho)
  outs = adjustCounts(sc)
  seurat <- CreateSeuratObject(counts = outs)
  # save
  save_path = file.path('/niuyw-usb-disk/Projects/scESCA/cellranger/soupx', sam)
  dir.create(save_path, showWarnings = FALSE)
  saveRDS(seurat, file = file.path(save_path, 'seurat.rds'))
}

# use mast markers to estimate soup
sams = c('S0619_N1', 'S1204_N5')
for (sam in sams) {
  data_dir = file.path(cellranger_dir, sam)
  raw_count_dir = file.path(data_dir, 'raw_feature_bc_matrix')
  filtered_count_dir = file.path(data_dir, 'filtered_feature_bc_matrix')
  # prepare for soupx
  seurat.raw = .build_seurat(filtered_count_dir, sam)
  seurat = .filter_cells(seurat.raw)
  seurat

  tmp = .soupx(raw_count_dir, seurat)
  sc = tmp$sc
  DR = tmp$DR
  #plotMarkerDistribution(sc)

  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = soupx_genes$mast))
  #plotMarkerMap(sc, geneSet = soupx_genes$mast, DR = DR, useToEst = useToEst)

  sc = calculateContaminationFraction(sc, list(IG = soupx_genes$mast), useToEst = useToEst)
  mean(sc$metaData$rho)
  outs = adjustCounts(sc)
  seurat <- CreateSeuratObject(counts = outs)
  # save
  save_path = file.path('/niuyw-usb-disk/Projects/scESCA/cellranger/soupx', sam)
  dir.create(save_path, showWarnings = FALSE)
  saveRDS(seurat, file = file.path(save_path, 'seurat.rds'))
}

# use epithelial markers to estimate soup
sams = c('S0619_N2', 'S0730_N1', 'S0730_N2', 'S0920_N1', 'S0920_N2', 'S1125_N1')
for (sam in sams) {
  data_dir = file.path(cellranger_dir, sam)
  raw_count_dir = file.path(data_dir, 'raw_feature_bc_matrix')
  filtered_count_dir = file.path(data_dir, 'filtered_feature_bc_matrix')
  # prepare for soupx
  seurat.raw = .build_seurat(filtered_count_dir, sam)
  seurat = .filter_cells(seurat.raw)
  seurat

  tmp = .soupx(raw_count_dir, seurat)
  sc = tmp$sc
  DR = tmp$DR
  #plotMarkerDistribution(sc)

  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = soupx_genes$epithelia))
  #plotMarkerMap(sc, geneSet = soupx_genes$epithelia, DR = DR, useToEst = useToEst)

  sc = calculateContaminationFraction(sc, list(IG = soupx_genes$epithelia), useToEst = useToEst)
  mean(sc$metaData$rho)
  outs = adjustCounts(sc)
  seurat <- CreateSeuratObject(counts = outs)
  # save
  save_path = file.path('/niuyw-usb-disk/Projects/scESCA/cellranger/soupx', sam)
  dir.create(save_path, showWarnings = FALSE)
  saveRDS(seurat, file = file.path(save_path, 'seurat.rds'))
}

# merge samples
sams = c('S0619_P1', 'S0619_P2', 'S0619_LN1', 'S0619_LN2', 'S0619_LN3', 'S0619_LN4', 'S0619_LN5', 'S0619_LN6', 'S0619_N1', 'S0619_N2', 'S0619_N3', 'S0619_N4', 'S0619_N5', 'S0730_P1', 'S0730_P2', 'S0730_LN1', 'S0730_LN2', 'S0730_LN3', 'S0730_LN4', 'S0730_N1', 'S0730_N2', 'S0730_N3', 'S0730_N4', 'S0730_N5', 'S0819_P1', 'S0819_P2', 'S0819_LN1', 'S0819_LN2', 'S0819_LN3', 'S0819_LN4', 'S0819_LN5', 'S0819_N1', 'S0819_N2', 'S0819_N3', 'S0819_N4', 'S0819_N5', 'S0920_P1', 'S0920_P2', 'S0920_LN1', 'S0920_LN2', 'S0920_LN3', 'S0920_LN4', 'S0920_N1', 'S0920_N2', 'S0920_N3', 'S0920_N4', 'S0920_N5', 'S1125_P1', 'S1125_P2', 'S1125_LN1', 'S1125_LN2', 'S1125_LN3', 'S1125_LN4', 'S1125_LN5', 'S1125_LN6', 'S1125_N1', 'S1125_N2', 'S1125_N3', 'S1125_N4', 'S1125_N5', 'S1204_P1', 'S1204_P2', 'S1204_LN1', 'S1204_LN2', 'S1204_LN3', 'S1204_LN5', 'S1204_LN6', 'S1204_LN7', 'S1204_N1', 'S1204_N2', 'S1204_N3', 'S1204_N4', 'S1204_N5')
datalist = list()
for (i in 1:length(sams)) {
  print(i)
  data_dir = file.path('/niuyw-usb-disk/Projects/scESCA/cellranger/soupx', sams[i])
  sam = str_replace(sams[i], '_', '.')
  seurat = readRDS(file.path(data_dir, 'seurat.rds'))
  seurat$orig.ident = sam
  seurat = RenameCells(seurat, new.names = str_replace(colnames(seurat), '$', paste0('-', i)))
  datalist[[sam]] = seurat
}

seurat = purrr::reduce(datalist, merge)

# save
DOCNAME = "01.preprocess_2"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)
#dir.create(here::here("docs", "assets", DOCNAME), showWarnings = FALSE)

saveRDS(seurat, file = here::here('output', DOCNAME, 'seurat.rds'))

