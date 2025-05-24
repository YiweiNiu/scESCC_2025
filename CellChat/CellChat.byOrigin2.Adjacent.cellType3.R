#!/usr/bin/env Rscript

# nohup Rscript analysis/CellChat.byOrigin2.Adjacent.cellType3.R > analysis/CellChat.byOrigin2.Adjacent.cellType3.log &

DOCNAME = "CellChat.byOrigin2.Adjacent.cellType3"
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
theme_set(theme_cowplot(font_size = 12,
                        rel_small = 10/12,
                        rel_tiny = 8/12,
                        rel_large = 12/12,
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

# cellchat
library(CellChat)
library(NMF)
library(ggalluvial)

# Seurat
library(Seurat)

options(stringsAsFactors = FALSE)

source(here::here("code/preprocess.R"))
source(here::here("code/glm.R"))
source(here::here("code/color.R"))
source(here::here("code/plot.R"))

seurat = readRDS(here::here('output/04.rm_cells/seurat.rds'))
seurat

srat_sub = subset(seurat, Origin2_n == "Adjacent")
rm(seurat)
data.input <- GetAssayData(srat_sub, assay = "RNA", slot = "data") # Tumorized data matrix
labels <- srat_sub$cellType3
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

# create
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# set the used database in the object
cellchat@DB <- CellChatDB

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

# Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")

# save
saveRDS(cellchat, file = here::here("output", DOCNAME, 'cellchat.rds'))

