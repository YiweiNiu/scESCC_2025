#!/usr/bin/env Rscript

# R 432
# nohup Rscript miloR.R > miloR.log &

word_dir <- "/work/home/project/scESCA/miloR_241215"
.libPaths('/work/home/project/20231127_DevM/devm_rproj/renv/library/R-4.3/x86_64-pc-linux-gnu')

DOCNAME <- "miloR"
dir.create(file.path(word_dir, DOCNAME), showWarnings = FALSE)


# Tidyverse
library(tidyverse)

library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)

# load
srat = readRDS('/work/home/project/scESCA/200227_6samples/output/04.rm_cells/seurat.rds')
srat = UpdateSeuratObject(srat)
srat$celltype <- srat$cellType3

# Create a Milo object
sce <- as.SingleCellExperiment(srat)
milo.obj <- Milo(sce)
milo.meta = srat@meta.data %>%
  rownames_to_column("barcode") %>%
  dplyr::select(barcode, Tissue = Origin2_n, Patient, Sample = Source) %>%
  mutate(across(everything(), as.character))

# Next we build the KNN graph and define neighbourhoods to quantify cell abundance across our experimental samples.
k = 15
prop = 0.05
use_rep = "PCA"
ncomps <- ncol(reducedDims(milo.obj)[[use_rep]])
milo.obj <- buildGraph(milo.obj, k = k, d = ncomps, reduced.dim = use_rep)
milo.obj <- makeNhoods(milo.obj, prop=prop, k=k, d=ncomps, reduced_dims = use_rep,
                       refined=TRUE, refinement_scheme="graph")

# size
median(colSums(nhoods(milo.obj)))
mean(colSums(nhoods(milo.obj)))
range(colSums(nhoods(milo.obj)))

# Counting cells in neighbourhoods
milo.obj <- countCells(milo.obj, samples="Sample", meta.data = milo.meta)

# Visualize neighbourhoods displaying DA
milo.obj <- buildNhoodGraph(milo.obj)


# We will use these contrasts to explicitly define which groups will be compared to each other.
thy.design <- milo.meta[,c("Tissue", "Patient", "Sample")]
thy.design <- distinct(thy.design)
rownames(thy.design) <- thy.design$Sample
## Reorder rownames to match columns of nhoodCounts(milo)
thy.design <- thy.design[colnames(nhoodCounts(milo.obj)), , drop=FALSE]
table(thy.design$Tissue)

# Make contrast
contrast.all <- c("TissuepostPBMC - TissueprePBMC",
                  "TissuepLN - TissuenLN",
                  "TissueAdjacent - TissueNormal",
                  "TissueTumor - TissueAdjacent",
                  "TissueTumor - TissueNormal")

# this is the edgeR code called by `testNhoods`
model <- model.matrix(~ 0 + Tissue, data=thy.design)
mod.constrast <- makeContrasts(contrasts=contrast.all, levels=model)
mod.constrast

# test
da_results <- testNhoods(milo.obj, design = ~ 0 + Tissue, design.df = thy.design,
                         model.contrasts = contrast.all,
                         fdr.weighting = "graph-overlap",
                         reduced.dim = use_rep)
table(da_results$SpatialFDR < 0.1)

# annotate
da_results <- annotateNhoods(milo.obj, da_results, coldata_col = "cellType3")

# save
saveRDS(da_results, file = file.path(word_dir, DOCNAME, "milo_res_byOrigin2_n.rds"))
saveRDS(milo.meta, file = file.path(word_dir, DOCNAME, "milo_meta.rds"))
saveRDS(thy.design, file = file.path(word_dir, DOCNAME, "milo_design.rds"))
saveRDS(milo.obj, file = file.path(word_dir, DOCNAME, "milo_obj.rds"))

Sys.Date()
sessionInfo()

