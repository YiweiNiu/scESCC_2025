#!/usr/bin/env Rscript

# source ~/software/anaconda3/bin/activate
# conda activate R-4.0.3
# nohup Rscript analysis/CytoTRACE.endothelia.R > analysis/CytoTRACE.endothelia.log &

DOCNAME = "CytoTRACE.endothelia"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

# Seurat
library(Seurat)

# CytoTRACE
library(CytoTRACE)

seurat = readRDS(here::here('output/04.rm_cells/seurat_endothelia.rds'))
Idents(seurat) <- 'seurat_clusters'
seurat

mat = GetAssayData(seurat, assay = 'RNA', slot = 'counts')
mat = Matrix::as.matrix(mat)
result = CytoTRACE(mat, enableFast = F, ncores = 8)
saveRDS(result, file = here::here("output", DOCNAME, "CytoTRACE.res.rds"))

# blood
seurat_subset = subset(seurat, subset = level_1 != "lymphatic")
mat = GetAssayData(seurat_subset, assay = 'RNA', slot = 'counts')
mat = Matrix::as.matrix(mat)
result = CytoTRACE(mat, enableFast = F, ncores = 10)
saveRDS(result, file = here::here("output", DOCNAME, "blood_endo.CytoTRACE.res.rds"))

