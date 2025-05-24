#!/usr/bin/env Rscript

# cd /niuyw-usb-disk/Projects/scESCA/200227_6samples
# nohup Rscript analysis/exp_similarity.cross_tissue.Origin2.R > analysis/exp_similarity.cross_tissue.Origin2.log &

DOCNAME = "exp_similarity.cross_tissue.Origin2"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

# Tidyverse
library(tidyverse)
# Plotting
library(ggplotify)
library(ggcorrplot)
library(ggpubr)
# patch
library(patchwork)
library(cowplot)
theme_set(theme_cowplot(font_size = 10,
                        rel_small = 8/10,
                        rel_tiny = 6/10,
                        rel_large = 10/10,
                        font_family = "Arial"))
# heatmap
library(pheatmap)
library(ComplexHeatmap)
# fonts
library(extrafont)
#font_import()
#font_import(paths = "/T02Data/niuyw/.local/share/fonts")
loadfonts(device = "pdf", quiet = TRUE)
# color
library(ggsci)
# BiocNeighbors
library(BiocNeighbors)
# Seurat
library(Seurat)

source(here::here("code/preprocess.R"))
source(here::here("code/glm.R"))
source(here::here("code/color.R"))
source(here::here("code/plot.R"))

seurat = readRDS(here::here('output/04.rm_cells/seurat.rds'))
Idents(seurat) <- 'cellType3'
seurat


# prePBMC
t = "prePBMC"
# Get matrix of 50 PCs.
srat_target <- subset(seurat, subset = Origin2_n == t)
srat_others <- subset(seurat, Origin2_n != t)
# get PCA score of each cell
target_matrix <- srat_target@reductions$pca[[]]
others_matrix <- srat_others@reductions$pca[[]]
KNN_result <- queryKNN(others_matrix, target_matrix, k=1, BNPARAM = KmknnParam())
# generate tissue origin data frame and its summary table
tissue_origin <- as.data.frame(KNN_result$index)
colnames(tissue_origin) = "NN_index"
tissue_origin$NN_tissue = srat_others$Origin2_n[tissue_origin$NN_index]
tissue_origin$NN_type = srat_others$cellType3[tissue_origin$NN_index]
# the cell type
rownames(tissue_origin) = rownames(target_matrix)
tissue_origin$cellType3 = srat_target$cellType3
# change order of columns
tissue_origin = tissue_origin[,c("cellType3", "NN_index", "NN_tissue", "NN_type")]
# test
set.seed(123)
fake_tissue_origin = sapply(1:1000, function(x){
  sample(srat_others$Origin2_n)[tissue_origin$NN_index]
})
fake_tissue_origin = as.data.frame(fake_tissue_origin)
fake_tissue_origin[1:3, 1:3]
# observed
real_num = table(tissue_origin$cellType3, as.character(tissue_origin$NN_tissue))
# simulated distribution
fake_num = lapply(fake_tissue_origin, FUN = function(x){
  table(tissue_origin$cellType3, x)
})
# count the number when the permutated num over the observed num
tmp_count = lapply(fake_num, FUN = function(x){
  a = x - real_num
  a[a < 0] = 0
  a[a > 0] = 1
  a
})
tot_count = Reduce("+", tmp_count)
# empirical p value
p_mat = tot_count/1000
# save
saveRDS(tissue_origin, file = here::here("output", DOCNAME, paste0(t, ".tissue_origin.rds")))
saveRDS(p_mat, file = here::here("output", DOCNAME, paste0(t, ".p_mat.rds")))

# postPBMC
t = "postPBMC"
# Get matrix of 50 PCs.
srat_target <- subset(seurat, subset = Origin2_n == t)
srat_others <- subset(seurat, Origin2_n != t)
# get PCA score of each cell
target_matrix <- srat_target@reductions$pca[[]]
others_matrix <- srat_others@reductions$pca[[]]
KNN_result <- queryKNN(others_matrix, target_matrix, k=1, BNPARAM = KmknnParam())
# generate tissue origin data frame and its summary table
tissue_origin <- as.data.frame(KNN_result$index)
colnames(tissue_origin) = "NN_index"
tissue_origin$NN_tissue = srat_others$Origin2_n[tissue_origin$NN_index]
tissue_origin$NN_type = srat_others$cellType3[tissue_origin$NN_index]
# the cell type
rownames(tissue_origin) = rownames(target_matrix)
tissue_origin$cellType3 = srat_target$cellType3
# change order of columns
tissue_origin = tissue_origin[,c("cellType3", "NN_index", "NN_tissue", "NN_type")]
# test
set.seed(123)
fake_tissue_origin = sapply(1:1000, function(x){
  sample(srat_others$Origin2_n)[tissue_origin$NN_index]
})
fake_tissue_origin = as.data.frame(fake_tissue_origin)
fake_tissue_origin[1:3, 1:3]
# observed
real_num = table(tissue_origin$cellType3, as.character(tissue_origin$NN_tissue))
# simulated distribution
fake_num = lapply(fake_tissue_origin, FUN = function(x){
  table(tissue_origin$cellType3, x)
})
# count the number when the permutated num over the observed num
tmp_count = lapply(fake_num, FUN = function(x){
  a = x - real_num
  a[a < 0] = 0
  a[a > 0] = 1
  a
})
tot_count = Reduce("+", tmp_count)
# empirical p value
p_mat = tot_count/1000
# save
saveRDS(tissue_origin, file = here::here("output", DOCNAME, paste0(t, ".tissue_origin.rds")))
saveRDS(p_mat, file = here::here("output", DOCNAME, paste0(t, ".p_mat.rds")))

# nLN
t = "nLN"
# Get matrix of 50 PCs.
srat_target <- subset(seurat, subset = Origin2_n == t)
srat_others <- subset(seurat, Origin2_n != t)
# get PCA score of each cell
target_matrix <- srat_target@reductions$pca[[]]
others_matrix <- srat_others@reductions$pca[[]]
KNN_result <- queryKNN(others_matrix, target_matrix, k=1, BNPARAM = KmknnParam())
# generate tissue origin data frame and its summary table
tissue_origin <- as.data.frame(KNN_result$index)
colnames(tissue_origin) = "NN_index"
tissue_origin$NN_tissue = srat_others$Origin2_n[tissue_origin$NN_index]
tissue_origin$NN_type = srat_others$cellType3[tissue_origin$NN_index]
# the cell type
rownames(tissue_origin) = rownames(target_matrix)
tissue_origin$cellType3 = srat_target$cellType3
# change order of columns
tissue_origin = tissue_origin[,c("cellType3", "NN_index", "NN_tissue", "NN_type")]
# test
set.seed(123)
fake_tissue_origin = sapply(1:1000, function(x){
  sample(srat_others$Origin2_n)[tissue_origin$NN_index]
})
fake_tissue_origin = as.data.frame(fake_tissue_origin)
fake_tissue_origin[1:3, 1:3]
# observed
real_num = table(tissue_origin$cellType3, as.character(tissue_origin$NN_tissue))
# simulated distribution
fake_num = lapply(fake_tissue_origin, FUN = function(x){
  table(tissue_origin$cellType3, x)
})
# count the number when the permutated num over the observed num
tmp_count = lapply(fake_num, FUN = function(x){
  a = x - real_num
  a[a < 0] = 0
  a[a > 0] = 1
  a
})
tot_count = Reduce("+", tmp_count)
# empirical p value
p_mat = tot_count/1000
# save
saveRDS(tissue_origin, file = here::here("output", DOCNAME, paste0(t, ".tissue_origin.rds")))
saveRDS(p_mat, file = here::here("output", DOCNAME, paste0(t, ".p_mat.rds")))


# pLN
t = "pLN"
# Get matrix of 50 PCs.
srat_target <- subset(seurat, subset = Origin2_n == t)
srat_others <- subset(seurat, Origin2_n != t)
# get PCA score of each cell
target_matrix <- srat_target@reductions$pca[[]]
others_matrix <- srat_others@reductions$pca[[]]
KNN_result <- queryKNN(others_matrix, target_matrix, k=1, BNPARAM = KmknnParam())
# generate tissue origin data frame and its summary table
tissue_origin <- as.data.frame(KNN_result$index)
colnames(tissue_origin) = "NN_index"
tissue_origin$NN_tissue = srat_others$Origin2_n[tissue_origin$NN_index]
tissue_origin$NN_type = srat_others$cellType3[tissue_origin$NN_index]
# the cell type
rownames(tissue_origin) = rownames(target_matrix)
tissue_origin$cellType3 = srat_target$cellType3
# change order of columns
tissue_origin = tissue_origin[,c("cellType3", "NN_index", "NN_tissue", "NN_type")]
# test
set.seed(123)
fake_tissue_origin = sapply(1:1000, function(x){
  sample(srat_others$Origin2_n)[tissue_origin$NN_index]
})
fake_tissue_origin = as.data.frame(fake_tissue_origin)
fake_tissue_origin[1:3, 1:3]
# observed
real_num = table(tissue_origin$cellType3, as.character(tissue_origin$NN_tissue))
# simulated distribution
fake_num = lapply(fake_tissue_origin, FUN = function(x){
  table(tissue_origin$cellType3, x)
})
# count the number when the permutated num over the observed num
tmp_count = lapply(fake_num, FUN = function(x){
  a = x - real_num
  a[a < 0] = 0
  a[a > 0] = 1
  a
})
tot_count = Reduce("+", tmp_count)
# empirical p value
p_mat = tot_count/1000
# save
saveRDS(tissue_origin, file = here::here("output", DOCNAME, paste0(t, ".tissue_origin.rds")))
saveRDS(p_mat, file = here::here("output", DOCNAME, paste0(t, ".p_mat.rds")))

# Normal
t = "Normal"
# Get matrix of 50 PCs.
srat_target <- subset(seurat, subset = Origin2_n == t)
srat_others <- subset(seurat, Origin2_n != t)
# get PCA score of each cell
target_matrix <- srat_target@reductions$pca[[]]
others_matrix <- srat_others@reductions$pca[[]]
KNN_result <- queryKNN(others_matrix, target_matrix, k=1, BNPARAM = KmknnParam())
# generate tissue origin data frame and its summary table
tissue_origin <- as.data.frame(KNN_result$index)
colnames(tissue_origin) = "NN_index"
tissue_origin$NN_tissue = srat_others$Origin2_n[tissue_origin$NN_index]
tissue_origin$NN_type = srat_others$cellType3[tissue_origin$NN_index]
# the cell type
rownames(tissue_origin) = rownames(target_matrix)
tissue_origin$cellType3 = srat_target$cellType3
# change order of columns
tissue_origin = tissue_origin[,c("cellType3", "NN_index", "NN_tissue", "NN_type")]
# test
set.seed(123)
fake_tissue_origin = sapply(1:1000, function(x){
  sample(srat_others$Origin2_n)[tissue_origin$NN_index]
})
fake_tissue_origin = as.data.frame(fake_tissue_origin)
fake_tissue_origin[1:3, 1:3]
# observed
real_num = table(tissue_origin$cellType3, as.character(tissue_origin$NN_tissue))
# simulated distribution
fake_num = lapply(fake_tissue_origin, FUN = function(x){
  table(tissue_origin$cellType3, x)
})
# count the number when the permutated num over the observed num
tmp_count = lapply(fake_num, FUN = function(x){
  a = x - real_num
  a[a < 0] = 0
  a[a > 0] = 1
  a
})
tot_count = Reduce("+", tmp_count)
# empirical p value
p_mat = tot_count/1000
# save
saveRDS(tissue_origin, file = here::here("output", DOCNAME, paste0(t, ".tissue_origin.rds")))
saveRDS(p_mat, file = here::here("output", DOCNAME, paste0(t, ".p_mat.rds")))

# Adjacent
t = "Adjacent"
# Get matrix of 50 PCs.
srat_target <- subset(seurat, subset = Origin2_n == t)
srat_others <- subset(seurat, Origin2_n != t)
# get PCA score of each cell
target_matrix <- srat_target@reductions$pca[[]]
others_matrix <- srat_others@reductions$pca[[]]
KNN_result <- queryKNN(others_matrix, target_matrix, k=1, BNPARAM = KmknnParam())
# generate tissue origin data frame and its summary table
tissue_origin <- as.data.frame(KNN_result$index)
colnames(tissue_origin) = "NN_index"
tissue_origin$NN_tissue = srat_others$Origin2_n[tissue_origin$NN_index]
tissue_origin$NN_type = srat_others$cellType3[tissue_origin$NN_index]
# the cell type
rownames(tissue_origin) = rownames(target_matrix)
tissue_origin$cellType3 = srat_target$cellType3
# change order of columns
tissue_origin = tissue_origin[,c("cellType3", "NN_index", "NN_tissue", "NN_type")]
# test
set.seed(123)
fake_tissue_origin = sapply(1:1000, function(x){
  sample(srat_others$Origin2_n)[tissue_origin$NN_index]
})
fake_tissue_origin = as.data.frame(fake_tissue_origin)
fake_tissue_origin[1:3, 1:3]
# observed
real_num = table(tissue_origin$cellType3, as.character(tissue_origin$NN_tissue))
# simulated distribution
fake_num = lapply(fake_tissue_origin, FUN = function(x){
  table(tissue_origin$cellType3, x)
})
# count the number when the permutated num over the observed num
tmp_count = lapply(fake_num, FUN = function(x){
  a = x - real_num
  a[a < 0] = 0
  a[a > 0] = 1
  a
})
tot_count = Reduce("+", tmp_count)
# empirical p value
p_mat = tot_count/1000
# save
saveRDS(tissue_origin, file = here::here("output", DOCNAME, paste0(t, ".tissue_origin.rds")))
saveRDS(p_mat, file = here::here("output", DOCNAME, paste0(t, ".p_mat.rds")))


# Tumor
t = "Tumor"
# Get matrix of 50 PCs.
srat_target <- subset(seurat, subset = Origin2_n == t)
srat_others <- subset(seurat, Origin2_n != t)
# get PCA score of each cell
target_matrix <- srat_target@reductions$pca[[]]
others_matrix <- srat_others@reductions$pca[[]]
KNN_result <- queryKNN(others_matrix, target_matrix, k=1, BNPARAM = KmknnParam())
# generate tissue origin data frame and its summary table
tissue_origin <- as.data.frame(KNN_result$index)
colnames(tissue_origin) = "NN_index"
tissue_origin$NN_tissue = srat_others$Origin2_n[tissue_origin$NN_index]
tissue_origin$NN_type = srat_others$cellType3[tissue_origin$NN_index]
# the cell type
rownames(tissue_origin) = rownames(target_matrix)
tissue_origin$cellType3 = srat_target$cellType3
# change order of columns
tissue_origin = tissue_origin[,c("cellType3", "NN_index", "NN_tissue", "NN_type")]
# test
set.seed(123)
fake_tissue_origin = sapply(1:1000, function(x){
  sample(srat_others$Origin2_n)[tissue_origin$NN_index]
})
fake_tissue_origin = as.data.frame(fake_tissue_origin)
fake_tissue_origin[1:3, 1:3]
# observed
real_num = table(tissue_origin$cellType3, as.character(tissue_origin$NN_tissue))
# simulated distribution
fake_num = lapply(fake_tissue_origin, FUN = function(x){
  table(tissue_origin$cellType3, x)
})
# count the number when the permutated num over the observed num
tmp_count = lapply(fake_num, FUN = function(x){
  a = x - real_num
  a[a < 0] = 0
  a[a > 0] = 1
  a
})
tot_count = Reduce("+", tmp_count)
# empirical p value
p_mat = tot_count/1000
# save
saveRDS(tissue_origin, file = here::here("output", DOCNAME, paste0(t, ".tissue_origin.rds")))
saveRDS(p_mat, file = here::here("output", DOCNAME, paste0(t, ".p_mat.rds")))

Sys.Date()
sessionInfo()

