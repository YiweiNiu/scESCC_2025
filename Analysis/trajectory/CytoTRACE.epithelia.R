#!/usr/bin/env Rscript

# nohup Rscript analysis/CytoTRACE.epithelia.R > analysis/CytoTRACE.epithelia.log &

DOCNAME = "CytoTRACE.epithelia"
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

seurat = readRDS(here::here('output/04.rm_cells/seurat_epithelia.rds'))
Idents(seurat) <- 'seurat_clusters'
seurat

mat = GetAssayData(seurat, assay = 'RNA', slot = 'counts')
mat = Matrix::as.matrix(mat)
result = CytoTRACE(mat, enableFast = F, ncores = 8)
saveRDS(result, file = here::here("output", DOCNAME, "CytoTRACE.res.rds"))
rm(seurat)

patients = c("S0619", "S0730", "S0819", "S0920", "S1125", "S1204")
lapply(patients, function(x){
  srat = readRDS(here::here(paste0(paste0('output/04.rm_cells/seurat_epithelia.', x), '.rds')))
  mat = GetAssayData(srat, assay = 'RNA', slot = 'counts')
  mat = Matrix::as.matrix(mat)
  res = CytoTRACE(mat, enableFast = F, ncores = 8)
  saveRDS(res, file = here::here("output", DOCNAME, paste0(paste0('CytoTRACE.res.', x), '.rds')))
})

