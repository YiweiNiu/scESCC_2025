#!/usr/bin/env Rscript

library(tidyverse)

# Seurat
library(Seurat)

library(future)
plan("multiprocess", workers = 4)
# 20*1024*1024*1024
options(future.globals.maxSize = 21474836480)

source(here::here("code/preprocess.R"))
source(here::here("code/plot.R"))

# 这个脚本用来删除细胞

DOCNAME = "04.rm_cells"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

# T, pre: 209384, post: 200025
seurat = readRDS(here::here('output/02.cell_type/seurat_tcells.rds'))
seurat
table(seurat$RNA_snn_res.2)
c_2_rm = as.character(c(40, 41, 45, 47, 48, 53))
seurat = subset(seurat, subset = RNA_snn_res.2 %in% c_2_rm, invert = TRUE)
seurat
table(seurat$RNA_snn_res.2)
seurat = get_cluster(seurat)
seurat = seurat$seurat_obj
saveRDS(seurat, here::here('output/04.rm_cells/seurat_tcells.rds'))

# B, pre: 57504, post: 57211
seurat = readRDS(here::here('output/02.cell_type/seurat_bcells.rds'))
seurat
table(seurat$RNA_snn_res.2)
c_2_rm = as.character(c(36, 40))
seurat = subset(seurat, subset = RNA_snn_res.2 %in% c_2_rm, invert = TRUE)
seurat
table(seurat$RNA_snn_res.2)
seurat = get_cluster(seurat)
seurat = seurat$seurat_obj
saveRDS(seurat, here::here('output/04.rm_cells/seurat_bcells.rds'))

# Fibroblasts, pre: 21352, post: 15686
seurat = readRDS(here::here('output/02.cell_type/seurat_fibroblasts.rds'))
seurat
table(seurat$RNA_snn_res.2)
c_2_rm = as.character(c(0, 1, 12, 16, 19, 26, 29))
seurat = subset(seurat, subset = RNA_snn_res.2 %in% c_2_rm, invert = TRUE)
seurat
table(seurat$RNA_snn_res.2)
seurat = get_cluster(seurat)
seurat = seurat$seurat_obj
saveRDS(seurat, here::here('output/04.rm_cells/seurat_fibroblasts.rds'))

# 200403, Fibroblasts, round 2, pre: 15686, post: 14439
seurat = readRDS(here::here('output/04.rm_cells/seurat_fibroblasts.rds'))
seurat
table(seurat$RNA_snn_res.2)
c_2_rm = as.character(c(18, 19, 22, 25, 27, 28, 29, 30))
seurat = subset(seurat, subset = RNA_snn_res.2 %in% c_2_rm, invert = TRUE)
seurat
table(seurat$RNA_snn_res.2)
seurat = get_cluster(seurat)
seurat = seurat$seurat_obj
saveRDS(seurat, here::here('output/04.rm_cells/seurat_fibroblasts.2.rds'))

# Endothelia, pre: 3010, post: 2580
seurat = readRDS(here::here('output/02.cell_type/seurat_endothelia.rds'))
seurat
table(seurat$RNA_snn_res.2)
c_2_rm = as.character(c(8, 15:19))
seurat = subset(seurat, subset = RNA_snn_res.2 %in% c_2_rm, invert = TRUE)
seurat
table(seurat$RNA_snn_res.2)
seurat = get_cluster(seurat)
seurat = seurat$seurat_obj
saveRDS(seurat, here::here('output/04.rm_cells/seurat_endothelia.rds'))

# Epithelia, pre: 61747, post: 37091
seurat = readRDS(here::here('output/02.cell_type/seurat_epithelia.rds'))
seurat
table(seurat$RNA_snn_res.2)
c_2_rm = as.character(c(2, 3, 6, 10, 11, 12, 15:17, 23, 25, 32, 33, 38:42, 45, 47:49))
seurat = subset(seurat, subset = RNA_snn_res.2 %in% c_2_rm, invert = TRUE)
seurat
table(seurat$RNA_snn_res.2)
seurat = get_cluster(seurat)
seurat = seurat$seurat_obj
saveRDS(seurat, here::here('output/04.rm_cells/seurat_epithelia.rds'))

# 200415, re-run Myeloid, pre: 25800, post: 21129
seurat = readRDS(here::here('output/02.cell_type/seurat_myeloid.rds'))
seurat
table(seurat$RNA_snn_res.2)
# cluster 40 was deleted due to low expression of lineage markers
# cluster 36 and cluster 43 were kept, because cluster 36 is a DC cluster that we want to keep
c_2_rm = as.character(c(7, 8, 19, 24, 27, 29, 35, 38, 39, 40, 42, 44:49))
seurat = subset(seurat, subset = RNA_snn_res.2 %in% c_2_rm, invert = TRUE)
seurat
table(seurat$RNA_snn_res.2)
seurat = get_cluster(seurat)
seurat = seurat$seurat_obj
saveRDS(seurat, here::here('output/04.rm_cells/seurat_myeloid.rds'))


