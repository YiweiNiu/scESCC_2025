#!/usr/bin/env Rscript

# source ~/software/anaconda3/bin/activate
# conda activate R-4.0.3
# nohup Rscript analysis/CytoTRACE.fibroblasts.R > analysis/CytoTRACE.fibroblasts.log &

DOCNAME = "CytoTRACE.fibroblasts"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)


# Seurat
library(Seurat)

# CytoTRACE
library(CytoTRACE)

seurat = readRDS(here::here('output/04.rm_cells/seurat_fibroblasts.rds'))
Idents(seurat) <- 'seurat_clusters'
seurat

mat = GetAssayData(seurat, assay = 'RNA', slot = 'counts')
mat = Matrix::as.matrix(mat)
result = CytoTRACE(mat, enableFast = F, ncores = 8)
saveRDS(result, file = here::here("output", DOCNAME, "CytoTRACE.res.rds"))

# pericyte + myCAF
seurat_subset = subset(seurat, subset = level_2 %in% c("myCAF", "pericyte"))
mat = GetAssayData(seurat_subset, assay = 'RNA', slot = 'counts')
mat = Matrix::as.matrix(mat)
result = CytoTRACE(mat, enableFast = F, ncores = 10)
saveRDS(result, file = here::here("output", DOCNAME, "myCAF_pericyte.CytoTRACE.res.rds"))

# NMF + iCAF
seurat_subset = subset(seurat, subset = level_2 %in% c("NMF", "iCAF"))
mat = GetAssayData(seurat_subset, assay = 'RNA', slot = 'counts')
mat = Matrix::as.matrix(mat)
result = CytoTRACE(mat, enableFast = F, ncores = 10)
saveRDS(result, file = here::here("output", DOCNAME, "NMF_iCAF.CytoTRACE.res.rds"))



