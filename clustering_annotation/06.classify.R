#!/usr/bin/env Rscript

# This script is to do the common tasks for cell classification

library(tidyverse)

# Plotting
library(cowplot)
library(pheatmap)
library(ggsci)
library(ggpubr)
library(ggcorrplot)
theme_set(theme_cowplot())

# Seurat
library(Seurat)
# sceasy
library(sceasy)

library(future)
plan("multiprocess", workers = 4)
# 20*1024*1024*1024
options(future.globals.maxSize = 21474836480)

source(here::here("code/preprocess.R"))


# T cells
seurat = readRDS(here::here('output/04.rm_cells/seurat_tcells.rds'))
Idents(seurat) <- 'seurat_clusters'
seurat

# level 1 sub-types: CD4/CD8/NK/Unknown
Idents(seurat) <- 'seurat_clusters'
seurat <- RenameIdents(seurat, `0` = "CD4", `1` = "CD4", `2` = "CD8", `3` = "CD4", `4` = "CD4", `5` = "CD4",
                       `6` = "CD8", `7` = "CD8", `8` = "CD4", `9` = "NK/NKT", `10` = "CD8",
                       `11` = "CD4", `12` = "Treg", `13` = "CD4", `14` = "CD8", `15` = "CD8",
                       `16` = "Treg", `17` = "CD4", `18` = "CD8", `19` = "NK/NKT", `20` = "CD4",
                       `21` = "CD8", `22` = "γδT", `23` = "CD8", `24` = "CD4", `25` = "CD8",
                       `26` = "Treg", `27` = "γδT", `28` = "CD8", `29` = "Unknown",
                       `30` = "CD8", `31` = "CD4")
# set level
Idents(seurat) = factor(Idents(seurat),
                        levels = c("CD4", "Treg", "CD8", "γδT", "NK/NKT", "Unknown"))
# add cell type to metadata
seurat$level_1 = Idents(seurat)

# level 2 sub-types: Naive/Tcm/Tem/Trm/Teff/Tex
Idents(seurat) <- 'seurat_clusters'
seurat <- RenameIdents(seurat, `0` = "Tcm", `1` = "Tn", `2` = "Tcm", `3` = "Tn", `4` = "Tn", `5` = "Tcm",
                       `6` = "Tex", `7` = "Tem", `8` = "Tn", `9` = "NK/NKT", `10` = "Teff",
                       `11` = "Tn", `12` = "Treg", `13` = "Tfh", `14` = "Teff", `15` = "Teff",
                       `16` = "Treg", `17` = "Tn", `18` = "Tex", `19` = "NK/NKT", `20` = "Tn",
                       `21` = "Tn", `22` = "γδT", `23` = "Trm", `24` = "Tem", `25` = "Teff",
                       `26` = "Treg", `27` = "γδT", `28` = "Tex", `29` = "Unknown", `30` = "Tex",
                       `31` = "Tcm")
# set level
Idents(seurat) = factor(Idents(seurat),
                        levels = c("Tn", "Tcm", "Tem", "Trm", "Tfh", "Teff", "Tex", "Treg", "γδT", "NK/NKT", "Unknown"))
# add cell type to metadata
seurat$level_2 = Idents(seurat)

# level 3 sub-types
Idents(seurat) <- 'seurat_clusters'
seurat <- RenameIdents(seurat, `0` = "CD4-C1-Tcm", `1` = "CD4-C2-Tn", `2` = "CD8-C1-Tcm", `3` = "CD4-C3-Tn",
                       `4` = "CD4-C4-Tn", `5` = "CD4-C5-Tcm", `6` = "CD8-C2-Tex", `7` = "CD8-C3-Tem",
                       `8` = "CD4-C6-Tn", `9` = "NK/NKT", `10` = "CD8-C4-Teff", `11` = "CD4-C7-Tn",
                       `12` = "Treg-C1", `13` = "CD4-C8-Tfh", `14` = "CD8-C5-Teff", `15` = "CD8-C6-Teff",
                       `16` = "Treg-C2", `17` = "CD4-C9-Tn", `18` = "CD8-C7-Tex", `19` = "NK/NKT",
                       `20` = "CD4-C10-Tn", `21` = "CD8-C8-Tn", `22` = "γδT-C1", `23` = "CD8-C9-Trm",
                       `24` = "CD4-C11-Tem", `25` = "CD8-C10-Teff", `26` = "Treg-C3",
                       `27` = "γδT-C2", `28` = "CD8-C11-Tex", `29` = "Unknown", `30` = "CD8-C12-Tex",
                       `31` = "CD4-C12-Tcm")
# set level
Idents(seurat) = factor(Idents(seurat),
                        levels = c("CD4-C1-Tcm", "CD4-C2-Tn", "CD4-C3-Tn", "CD4-C4-Tn", "CD4-C5-Tcm",
                                   "CD4-C6-Tn", "CD4-C7-Tn", "CD4-C8-Tfh", "CD4-C9-Tn", "CD4-C10-Tn",
                                   "CD4-C11-Tem", "CD4-C12-Tcm",
                                   # Treg
                                   "Treg-C1", "Treg-C2", "Treg-C3",
                                   # CD8
                                   "CD8-C1-Tcm", "CD8-C2-Tex", "CD8-C3-Tem", "CD8-C4-Teff", "CD8-C5-Teff",
                                   "CD8-C6-Teff", "CD8-C7-Tex", "CD8-C8-Tn", "CD8-C9-Trm", "CD8-C10-Teff",
                                   "CD8-C11-Tex", "CD8-C12-Tex",
                                   # γδT
                                   "γδT-C1", "γδT-C2",
                                   "NK/NKT", "Unknown"))
# add cell type to metadata
seurat$level_3 = Idents(seurat)

# save
Idents(seurat) <- 'seurat_clusters'
saveRDS(seurat, file = here::here('output/04.rm_cells/seurat_tcells.rds'))
# to h5ad
sceasy::convertFormat(seurat, from="seurat", to="anndata",
                      main_layer = 'counts',
                      outFile=here::here('output/04.rm_cells', paste0('seurat_tcells', '.h5ad')))
# export cellmeta
vars = c("UMAP_1", "UMAP_2", "seurat_clusters", 'cellType', 'Patient', 'Tissue',
         'Origin2_n', "Origin3", 'Metastasis_n', 'Source',
         "Origin2", "Origin", "Metastasis",
         "level_1", "level_2", "level_3")
metadata <- FetchData(seurat,
                      vars = vars) %>%
  rownames_to_column('barcode')
write_csv(metadata, file = 'output/04.rm_cells/seurat_tcells.cellmeta.csv')
# CD8+ T cells
seurat_cd8 = subset(seurat, subset = level_1 == "CD8")
seurat_cd8 = get_projection(seurat_cd8)
saveRDS(seurat_cd8, file = here::here('output/04.rm_cells/seurat_cd8.rds'))
# to h5ad
sceasy::convertFormat(seurat_cd8, from="seurat", to="anndata",
                      main_layer = 'counts',
                      outFile=here::here('output/04.rm_cells', paste0('seurat_cd8', '.h5ad')))
# CD4+ T cells
seurat_cd4 = subset(seurat, subset = level_1 == "CD4")
seurat_cd4 = get_projection(seurat_cd4)
saveRDS(seurat_cd4, file = here::here('output/04.rm_cells/seurat_cd4.rds'))
# to h5ad
sceasy::convertFormat(seurat_cd4, from="seurat", to="anndata",
                      main_layer = 'counts',
                      outFile=here::here('output/04.rm_cells', paste0('seurat_cd4', '.h5ad')))
# Treg
seurat_treg = subset(seurat, subset = level_1 == "Treg")
seurat_treg = get_projection(seurat_treg)
saveRDS(seurat_treg, file = here::here('output/04.rm_cells/seurat_treg.rds'))
sceasy::convertFormat(seurat_treg, from="seurat", to="anndata",
                      main_layer = 'counts',
                      outFile=here::here('output/04.rm_cells', paste0('seurat_treg', '.h5ad')))
# CD4 + Treg
seurat_cd4_treg = subset(seurat, subset = level_1 %in% c("CD4", "Treg"))
seurat_cd4_treg = get_projection(seurat_cd4_treg)
saveRDS(seurat_cd4_treg, file = here::here('output/04.rm_cells/seurat_cd4_treg.rds'))
vars = c("UMAP_1", "UMAP_2", "seurat_clusters", 'Tissue', 'Origin',
         "Origin2_n", 'Origin3', 'Origin4', "Metastatic",
         'Patient', 'Source', 'cellType', "level_1", "level_2", "level_3")
FetchData(seurat_cd4_treg, vars = vars) %>%
  rownames_to_column('barcode') %>%
  write_csv(file = 'output/04.rm_cells/seurat_cd4_treg.cellmeta.csv')

# gdT
seurat_gdt = subset(seurat, subset = level_1 == "γδT")
seurat_gdt = get_projection(seurat_gdt)
saveRDS(seurat_gdt, file = here::here('output/04.rm_cells/seurat_gdt.rds'))

# Myeloid
seurat = readRDS(here::here('output/04.rm_cells/seurat_myeloid.rds'))
Idents(seurat) <- 'seurat_clusters'
seurat

# level 1 sub-types: CD4/CD8/NK/Unknown
Idents(seurat) <- 'seurat_clusters'
seurat <- RenameIdents(seurat, `0` = "Mono", `1` = "MΦ", `2` = "Mono", `3` = "Mono", `4` = "MΦ",
                       `5` = "Mono", `6` = "Mono", `7` = "MΦ", `8` = "MΦ", `9` = "DC", `10` = "Mast",
                       `11` = "DC", `12` = "Mast", `13` = "DC", `14` = "Mast", `15` = "DC",
                       `16` = "MΦ", `17` = "MΦ", `18` = "DC", `19` = "Mono", `20` = "DC",
                       `21` = "DC", `22` = "MΦ", `23` = "CD34+ cells", `24` = "DC")
# set level
Idents(seurat) = factor(Idents(seurat),
                        levels = c("Mono", "MΦ", "DC", "Mast", "CD34+ cells"))
# add cell type to metadata
seurat$level_1 = Idents(seurat)

# level 2 sub-types: Naive/Tcm/Tem/Trm/Teff/Tex
Idents(seurat) <- 'seurat_clusters'
seurat <- RenameIdents(seurat, `0` = "cMono", `1` = "MΦ", `2` = "cMono", `3` = "cMono", `4` = "MΦ",
                       `5` = "cMono", `6` = "ncMono", `7` = "MΦ", `8` = "MΦ", `9` = "cDC2",
                       `10` = "Mast", `11` = "cDC2", `12` = "Mast", `13` = "tDC", `14` = "Mast",
                       `15` = "pDC", `16` = "MΦ", `17` = "MΦ", `18` = "cDC1", `19` = "cMono",
                       `20` = "cDC1", `21` = "pDC", `22` = "MΦ", `23` = "CD34+ cells", `24` = "pDC")
# set level
Idents(seurat) = factor(Idents(seurat),
                        levels = c("cMono", "ncMono", "MΦ", "cDC1", "cDC2", "tDC", "pDC", "Mast", "CD34+ cells"))
# add cell type to metadata
seurat$level_2 = Idents(seurat)

# level 3 sub-types
Idents(seurat) <- 'seurat_clusters'
seurat <- RenameIdents(seurat, `0` = "cMono-C1", `1` = "MΦ-C1", `2` = "cMono-C2", `3` = "cMono-C3",
                       `4` = "MΦ-C2", `5` = "cMono-C4", `6` = "ncMono-C5", `7` = "MΦ-C3",
                       `8` = "MΦ-C4", `9` = "cDC2-C1", `10` = "Mast-C1", `11` = "cDC2-C2",
                       `12` = "Mast-C2", `13` = "tDC-C3", `14` = "Mast-C3", `15` = "pDC-C4",
                       `16` = "MΦ-C5", `17` = "MΦ-C6", `18` = "cDC1-C5", `19` = "cMono-C6",
                       `20` = "cDC1-C6", `21` = "pDC-C7", `22` = "MΦ-C7", `23` = "CD34+ cells",
                       `24` = "pDC-C8")
# set level
Idents(seurat) = factor(Idents(seurat),
                        levels = c("cMono-C1", "cMono-C2", "cMono-C3", "cMono-C4", "ncMono-C5", "cMono-C6",
                                   "MΦ-C1", "MΦ-C2", "MΦ-C3", "MΦ-C4", "MΦ-C5", "MΦ-C6", "MΦ-C7",
                                   "cDC2-C1", "cDC2-C2", "tDC-C3", "pDC-C4", "cDC1-C5", "cDC1-C6", "pDC-C7", "pDC-C8",
                                   "Mast-C1", "Mast-C2", "Mast-C3",
                                   "CD34+ cells"))
# add cell type to metadata
seurat$level_3 = Idents(seurat)

# save
Idents(seurat) <- 'seurat_clusters'
saveRDS(seurat, file = here::here('output/04.rm_cells/seurat_myeloid.rds'))
# export cellmeta
vars = c("UMAP_1", "UMAP_2", "seurat_clusters", 'cellType', 'Patient', 'Tissue',
         'Origin2_n', "Origin3", 'Metastasis_n', 'Source',
         "Origin2", "Metastasis",
         "level_1", "level_2", "level_3")
metadata <- FetchData(seurat,
                      vars = vars) %>%
  rownames_to_column('barcode')
write_csv(metadata, file = 'output/04.rm_cells/seurat_myeloid.cellmeta.csv')
# Mono
seurat_mono = subset(seurat, subset = level_1 == "Mono")
seurat_mono = get_projection(seurat_mono)
saveRDS(seurat_mono, file = here::here('output/04.rm_cells/seurat_mono.rds'))
# MΦ
seurat_mac = subset(seurat, subset = level_1 == "MΦ")
seurat_mac = get_projection(seurat_mac)
saveRDS(seurat_mac, file = here::here('output/04.rm_cells/seurat_mac.rds'))
# Mono + MΦ
seurat_mps = subset(seurat, subset = level_1 %in% c("Mono", "MΦ"))
seurat_mps = get_projection(seurat_mps)
saveRDS(seurat_mps, file = here::here('output/04.rm_cells/seurat_mps.rds'))
# DC
seurat_dc = subset(seurat, subset = level_1 == "DC")
seurat_dc = get_projection(seurat_dc)
saveRDS(seurat_dc, file = here::here('output/04.rm_cells/seurat_dc.rds'))
# pDC
seurat_pDC = subset(seurat, subset = level_2 == "pDC")
seurat_pDC = get_projection(seurat_pDC)
saveRDS(seurat_pDC, file = here::here('output/04.rm_cells/seurat_pDC.rds'))
# mDC (cDC + tDC)
seurat_mDC = subset(seurat, subset = level_2 %in% c("cDC1", "cDC2", "tDC"))
seurat_mDC = get_projection(seurat_mDC)
saveRDS(seurat_mDC, file = here::here('output/04.rm_cells/seurat_mDC.rds'))
# cDC2 + tDC (20230621)
seurat_cDC2_tDC = subset(seurat, subset = level_2 %in% c("cDC2", "tDC"))
seurat_cDC2_tDC = get_projection(seurat_cDC2_tDC)
saveRDS(seurat_cDC2_tDC, file = here::here('output/04.rm_cells/seurat_cDC2_tDC.rds'))
# mast
seurat_mast = subset(seurat, subset = level_1 == "Mast")
seurat_mast = get_projection(seurat_mast)
saveRDS(seurat_mast, file = here::here('output/04.rm_cells/seurat_mast.rds'))


# Endothelia
seurat = readRDS(here::here('output/04.rm_cells/seurat_endothelia.rds'))
Idents(seurat) <- 'seurat_clusters'
seurat

# level 1 sub-types: blood/lymphatic
Idents(seurat) <- 'seurat_clusters'
seurat <- RenameIdents(seurat, `0` = "blood", `1` = "blood", `2` = "blood",
                       `3` = "blood", `4` = "lymphatic", `5` = "blood",
                       `6` = "lymphatic", `7` = "blood", `8` = "blood", `9` = "blood")
# set level
Idents(seurat) = factor(Idents(seurat),
                        levels = c("blood", "lymphatic"))
# add cell type to metadata
seurat$level_1 = Idents(seurat)

# level_2
Idents(seurat) <- 'seurat_clusters'
seurat <- RenameIdents(seurat, `0` = "immature", `1` = "pcv-C1", `2` = "pcv-C2",
                       `3` = "tip cell", `4` = "normal lymphatics",
                       `5` = "pcv-C3", `6` = "tumor lymphatics",
                       `7` = "arteries", `8` = "pcv-C4", `9` = "pcv-C5")
# set level
Idents(seurat) = factor(Idents(seurat),
                        levels = c("arteries", "pcv-C1", "pcv-C2", "pcv-C3",
                                   "pcv-C4", "pcv-C5", "immature", "tip cell",
                                   "normal lymphatics", "tumor lymphatics"))
# add cell type to metadata
seurat$level_2 = Idents(seurat)
# export cellmeta
vars = c("UMAP_1", "UMAP_2", "seurat_clusters", 'cellType', 'Patient', 'Tissue',
         'Origin2_n', "Origin3", 'Metastasis_n', 'Source',
         "Origin2", "Metastasis",
         "level_1", "level_2")
metadata <- FetchData(seurat,
                      vars = vars) %>%
  rownames_to_column('barcode')
write_csv(metadata, file = 'output/04.rm_cells/seurat_endothelia.cellmeta.csv')
saveRDS(seurat, file = here::here('output/04.rm_cells/seurat_endothelia.rds'))


# Fibroblasts
seurat = readRDS(here::here('output/04.rm_cells/seurat_fibroblasts.rds'))
seurat

# level_1 sub-types: fibroblasts/pericyte/VSMC
Idents(seurat) <- 'seurat_clusters'
seurat <- RenameIdents(seurat, `0` = "fibroblasts", `1` = "fibroblasts", `2` = "fibroblasts",
                       `3` = "fibroblasts", `4` = "pericyte", `5` = "fibroblasts",
                       `6` = "fibroblasts", `7` = "fibroblasts", `8` = "VSMC",
                       `9` = "fibroblasts", `10` = "pericyte", `11` = "fibroblasts",
                       `12` = "pericyte", `13` = "fibroblasts")
# set level
Idents(seurat) = factor(Idents(seurat),
                        levels = c("fibroblasts", "pericyte", "VSMC"))
# add cell type to metadata
seurat$level_1 = Idents(seurat)

# level_2 sub-types: apCAF/myCAF/iCAF
Idents(seurat) <- 'seurat_clusters'
seurat <- RenameIdents(seurat, `0` = "iCAF", `1` = "NMF", `2` = "myCAF",
                       `3` = "NMF", `4` = "pericyte", `5` = "myCAF",
                       `6` = "myCAF", `7` = "iCAF", `8` = "VSMC",
                       `9` = "iCAF", `10` = "pericyte", `11` = "NMF",
                       `12` = "pericyte", `13` = "apCAF")
# set level
Idents(seurat) = factor(Idents(seurat),
                        levels = c("NMF", "iCAF", "myCAF", "pericyte", "apCAF", "VSMC"))
# add cell type to metadata
seurat$level_2 = Idents(seurat)

# level_3 sub-types: cluster
Idents(seurat) <- 'seurat_clusters'
seurat <- RenameIdents(seurat, `0` = "iCAF-C1", `1` = "NMF-C1", `2` = "myCAF-C1",
                       `3` = "NMF-C2", `4` = "pericyte-C1", `5` = "myCAF-C2",
                       `6` = "myCAF-C3", `7` = "iCAF-C2", `8` = "VSMC",
                       `9` = "iCAF-C3", `10` = "pericyte-C2", `11` = "NMF-C3",
                       `12` = "pericyte-C3", `13` = "apCAF")
# set level
Idents(seurat) = factor(Idents(seurat),
                        levels = c("NMF-C1", "NMF-C2", "NMF-C3",
                                   "iCAF-C1", "iCAF-C2", "iCAF-C3",
                                   "myCAF-C1", "myCAF-C2", "myCAF-C3",
                                   "pericyte-C1", "pericyte-C2", "pericyte-C3",
                                   "apCAF",
                                   "VSMC"))
# add cell type to metadata
seurat$level_3 = Idents(seurat)
Idents(seurat) <- 'seurat_clusters'
# export cellmeta
vars = c("UMAP_1", "UMAP_2", "seurat_clusters", 'cellType', 'Patient', 'Tissue',
         'Origin2_n', "Origin3", 'Metastasis_n', 'Source',
         "Origin2", "Metastasis",
         "level_1", "level_2", "level_3")
metadata <- FetchData(seurat,
                      vars = vars) %>%
  rownames_to_column('barcode')
write_csv(metadata, path = 'output/04.rm_cells/seurat_fibroblasts.cellmeta.csv')
saveRDS(seurat, file = here::here('output/04.rm_cells/seurat_fibroblasts.rds'))


# B cells
seurat = readRDS(here::here('output/04.rm_cells/seurat_bcells.rds'))
Idents(seurat) <- 'seurat_clusters'
seurat

# level 1 sub-types: Naive/Mem/GCB/Plasma/DN/CD5
Idents(seurat) <- 'seurat_clusters'
seurat <- RenameIdents(seurat, `0` = "Mem", `1` = "Naive", `2` = "Mem", `3` = "Mem", `4` = "Mem",
                       `5` = "Naive", `6` = "PC", `7` = "Mem", `8` = "Naive", `9` = "CD5",
                       `10` = "GCB", `11` = "Naive", `12` = "PC", `13` = "Naive", `14` = "Naive",
                       `15` = "GCB", `16` = "CD5", `17` = "Naive", `18` = "Mem", `19` = "Naive",
                       `20` = "Mem", `21` = "Mem", `22` = "GCB", `23` = "DN", `24` = "DN",
                       `25` = "PC", `26` = "PC")
# set level
Idents(seurat) = factor(Idents(seurat),
                        levels = c("Naive", "Mem", "PC", "GCB", "DN", "CD5"))
# add cell type to metadata
seurat$level_1 = Idents(seurat)

# level_2 sub-types: cluster
Idents(seurat) <- 'seurat_clusters'
seurat <- RenameIdents(seurat, `0` = "Mem-C1", `1` = "Naive-C1", `2` = "Mem-C2", `3` = "Mem-C3",
                       `4` = "Mem-C4", `5` = "Naive-C2", `6` = "PC-C1", `7` = "Mem-C5", `8` = "Naive-C3",
                       `9` = "CD5-C1", `10` = "GCB-C1", `11` = "Naive-C4", `12` = "PC-C2", `13` = "Naive-C5",
                       `14` = "Naive-C6", `15` = "GCB-C2", `16` = "CD5-C2", `17` = "Naive-C7", `18` = "Mem-C6",
                       `19` = "Naive-C8", `20` = "Mem-C7", `21` = "Mem-C8", `22` = "GCB-C3", `23` = "DN-C1",
                       `24` = "DN-C2", `25` = "PC-C3", `26` = "PC-C4")
# set level
Idents(seurat) = factor(Idents(seurat),
                        levels = c("Naive-C1", "Naive-C2", "Naive-C3", "Naive-C4", "Naive-C5", "Naive-C6",
                                   "Naive-C7", "Naive-C8", "Mem-C1", "Mem-C2", "Mem-C3", "Mem-C4", "Mem-C5",
                                   "Mem-C6", "Mem-C7", "Mem-C8", "PC-C1", "PC-C2", "PC-C3", "PC-C4", "GCB-C1",
                                   "GCB-C2", "GCB-C3", "DN-C1", "DN-C2", "CD5-C1", "CD5-C2"))
# add cell type to metadata
seurat$level_2 = Idents(seurat)
# save
Idents(seurat) <- 'seurat_clusters'
saveRDS(seurat, file = here::here('output/04.rm_cells/seurat_bcells.rds'))
# to h5ad
sceasy::convertFormat(seurat, from="seurat", to="anndata",
                      main_layer = 'counts',
                      outFile=here::here('output/04.rm_cells', paste0('seurat_bcells', '.h5ad')))
# export cellmeta
vars = c("UMAP_1", "UMAP_2", "seurat_clusters", 'cellType', 'Patient', 'Tissue', "Origin",
         'Origin2_n', "Origin3", 'Metastasis_n', 'Source',
         "Origin2", "Metastasis",
         "level_1", "level_2")
metadata <- FetchData(seurat,
                      vars = vars) %>%
  rownames_to_column('barcode')
write_csv(metadata, file = 'output/04.rm_cells/seurat_bcells.cellmeta.csv')


# All
seurat = readRDS(here::here('output/04.rm_cells/seurat.rds'))
seurat
# cellType2, cellType3
meta.t = read_csv('output/04.rm_cells/seurat_tcells.cellmeta.csv') %>%
  dplyr::select(barcode,
                cellType2 = level_1,
                cellType3 = level_3)
meta.b = read_csv('output/04.rm_cells/seurat_bcells.cellmeta.csv') %>%
  dplyr::select(barcode,
                cellType2 = level_1,
                cellType3 = level_2) %>%
  mutate(cellType2 = paste0("B-", cellType2),
         cellType3 = paste0("B-", cellType3))
meta.mye = read_csv('output/04.rm_cells/seurat_myeloid.cellmeta.csv') %>%
  dplyr::select(barcode,
                cellType2 = level_1,
                cellType3 = level_3)
meta.endo = read_csv('output/04.rm_cells/seurat_endothelia.cellmeta.csv') %>%
  dplyr::select(barcode,
                cellType2 = level_1,
                cellType3 = level_2) %>%
  mutate(cellType2 = case_when(
    cellType2 == "blood" ~ "blood_endo",
    cellType2 == "lymphatic" ~ "lymphatic_endo"
  ))
meta.fib = read_csv('output/04.rm_cells/seurat_fibroblasts.cellmeta.csv') %>%
  dplyr::select(barcode,
                cellType2 = level_2,
                cellType3 = level_3)
meta.platelets = seurat@meta.data %>%
  rownames_to_column("barcode") %>%
  filter(cellType == 'Platelets') %>%
  dplyr::select(barcode,
                cellType2 = cellType,
                cellType3 = cellType)
meta.epi = read_csv('output/04.rm_cells/seurat_epithelia.cellmeta.csv') %>%
  dplyr::select(barcode,
                cellType2 = malignant_1,
                cellType3 = seurat_clusters) %>%
  mutate(cellType3 = paste0('Epi_', cellType3))
meta.all = do.call(rbind, list(meta.t, meta.b, meta.mye,
                               meta.epi, meta.fib, meta.endo,
                               meta.platelets))
seurat@meta.data = seurat@meta.data %>%
  rownames_to_column("barcode") %>%
  left_join(meta.all, by = 'barcode') %>%
  column_to_rownames("barcode")
table(seurat$cellType2, useNA = 'ifany')
table(seurat$cellType3, useNA = 'ifany')

# rds
saveRDS(seurat, file = here::here('output/04.rm_cells/seurat.rds'))

# metadata
vars = c("UMAP_1", "UMAP_2", "seurat_clusters", 'Patient', 'Tissue',
         'Origin2_n', "Origin3", 'Metastasis_n', 'Source',
         "Origin2", "Metastasis",
         'cellType', "cellType2", "cellType3")
metadata <- FetchData(seurat,
                      vars = vars) %>%
  rownames_to_column('barcode')
write_csv(metadata, file = 'output/04.rm_cells/seurat.cellmeta.csv')

# sam_info
write_csv(seurat@misc$sam_info, file = 'data/sample_info.2.csv')


#####################################################
# separate by Patient
#####################################################

patients = c("S0619", "S0730", "S0819", "S0920", "S1125", "S1204")
get_vars <- function(seurat_obj) {
  vst_each = c()
  source = FetchData(seurat_obj, vars = c('Source'))
  for (i in unique(seurat_obj$Source)) {
    x = seurat_obj[,which(source == i)]
    if (dim(x)[2] < 100) {
      next
    } else {
      x = FindVariableFeatures(x, selection.method = 'vst', nfeatures = 500)
      vst_each = c(VariableFeatures(x), vst_each)
    }
  }
  VariableFeatures(seurat_obj) = unique(vst_each)
  return(seurat_obj)
}
#p = 'S0619'
for (p in patients) {
  cat(p)
  # fetch one patient
  s = FetchData(object = seurat, vars = 'Patient')
  srat_subset <- seurat[, which(x = s$Patient == p)]
  # get pc
  srat_subset = get_vars(srat_subset)
  vars.to.regress = c('nCount_RNA', 'percent.mt', 'percent.ribo', 'percent.dissociation',
                      'percent.heat', 'CC.Difference')
  srat_subset <- ScaleData(srat_subset, features = VariableFeatures(srat_subset),
                           vars.to.regress = vars.to.regress, verbose = F)
  srat_subset = RunPCA(srat_subset, features = VariableFeatures(srat_subset), verbose = F)
  n_pc = min(determine_pc_num(srat_subset))
  # non-linear dimensional reduction (UMAP/tSNE)
  srat_subset = RunUMAP(srat_subset, dims = 1:n_pc, verbose = F)
  srat_subset = RunTSNE(srat_subset, dims = 1:n_pc, verbose = F)
  fname = paste0(paste0('seurat.', p), '.rds')
  saveRDS(srat_subset, file = here::here('output', '04.rm_cells', fname))
}


