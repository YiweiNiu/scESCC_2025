#!/usr/bin/env Rscript

# nohup Rscript analysis/CytoTRACE.cd8.R > analysis/CytoTRACE.cd8.log &

DOCNAME = "CytoTRACE.cd8"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

# Tidyverse
library(tidyverse)

# Plotting
library(ggplotify)
library(ggcorrplot)
library(cowplot)
library(patchwork)
library(ggpubr)
theme_set(theme_cowplot())

# heatmap
library(pheatmap)
library(ComplexHeatmap)

# color
library(ggsci)

# Seurat
library(Seurat)

# CytoTRACE
library(CytoTRACE)

# random1
seurat = readRDS(here::here('output/04.rm_cells/seurat_cd8.random1.rds'))
Idents(seurat) <- 'level_3'
mat = GetAssayData(seurat, assay = 'RNA', slot = 'counts')
mat = Matrix::as.matrix(mat)
result = CytoTRACE(mat, enableFast = F, ncores = 8)
saveRDS(result, file = here::here("output", DOCNAME, "CytoTRACE.random1.rds"))

# random2
seurat = readRDS(here::here('output/04.rm_cells/seurat_cd8.random2.rds'))
Idents(seurat) <- 'level_3'
mat = GetAssayData(seurat, assay = 'RNA', slot = 'counts')
mat = Matrix::as.matrix(mat)
result = CytoTRACE(mat, enableFast = F, ncores = 8)
saveRDS(result, file = here::here("output", DOCNAME, "CytoTRACE.random2.rds"))

# random3
seurat = readRDS(here::here('output/04.rm_cells/seurat_cd8.random3.rds'))
Idents(seurat) <- 'level_3'
mat = GetAssayData(seurat, assay = 'RNA', slot = 'counts')
mat = Matrix::as.matrix(mat)
result = CytoTRACE(mat, enableFast = F, ncores = 8)
saveRDS(result, file = here::here("output", DOCNAME, "CytoTRACE.random3.rds"))

# random4
seurat = readRDS(here::here('output/04.rm_cells/seurat_cd8.random4.rds'))
Idents(seurat) <- 'level_3'
mat = GetAssayData(seurat, assay = 'RNA', slot = 'counts')
mat = Matrix::as.matrix(mat)
result = CytoTRACE(mat, enableFast = F, ncores = 8)
saveRDS(result, file = here::here("output", DOCNAME, "CytoTRACE.random4.rds"))

# random5
seurat = readRDS(here::here('output/04.rm_cells/seurat_cd8.random5.rds'))
Idents(seurat) <- 'level_3'
mat = GetAssayData(seurat, assay = 'RNA', slot = 'counts')
mat = Matrix::as.matrix(mat)
result = CytoTRACE(mat, enableFast = F, ncores = 8)
saveRDS(result, file = here::here("output", DOCNAME, "CytoTRACE.random5.rds"))

# random6
seurat = readRDS(here::here('output/04.rm_cells/seurat_cd8.random6.rds'))
Idents(seurat) <- 'level_3'
mat = GetAssayData(seurat, assay = 'RNA', slot = 'counts')
mat = Matrix::as.matrix(mat)
result = CytoTRACE(mat, enableFast = F, ncores = 8)
saveRDS(result, file = here::here("output", DOCNAME, "CytoTRACE.random6.rds"))


