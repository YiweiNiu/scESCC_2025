#!/usr/bin/env Rscript

# R 432
# nohup Rscript miloR_optimize_k_p.R > miloR_optimize_k_p.log &

word_dir <- "/work/home/project/scESCA/miloR_241215"
.libPaths('/work/home/project/20231127_DevM/devm_rproj/renv/library/R-4.3/x86_64-pc-linux-gnu')

DOCNAME <- "miloR_optimize_k_p"
dir.create(file.path(word_dir, DOCNAME), showWarnings = FALSE)

library(tidyverse)

# plot
library(scattermore)
library(ggrastr)
library(ggbeeswarm)

remove_x_axis <- function() {
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
}

remove_y_axis <- function() {
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
}

# patch
library(patchwork)
library(cowplot)
theme_set(theme_cowplot(
  font_size = 10,
  rel_small = 8 / 10,
  rel_tiny = 8 / 10,
  rel_large = 10 / 10
))
mytheme <- theme_cowplot(
  font_size = 10,
  rel_small = 8 / 10,
  rel_tiny = 8 / 10,
  rel_large = 10 / 10
)

# color
library(ggsci)

# srat
library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)

# read
srat = readRDS('/work/home/project/scESCA/200227_6samples/output/04.rm_cells/seurat.rds')
srat = UpdateSeuratObject(srat)

# Create a Milo object
sce <- as.SingleCellExperiment(srat)
milo.obj <- Milo(sce)
milo.obj

build_grap <- function(milo_obj = NULL, k = NULL, use_rep = "PCA") {
  ncomps <- ncol(reducedDims(milo_obj)[[use_rep]])
  # Next we build the KNN graph and define neighbourhoods to quantify cell abundance across our experimental samples.
  milo_obj <- buildGraph(milo_obj, k = k, d = ncomps, reduced.dim = use_rep)
  return(milo_obj)
}

optimize_k_p <- function(milo_obj = NULL, k = NULL, p = NULL, use_rep = "PCA") {
  ncomps <- ncol(reducedDims(milo_obj)[[use_rep]])
  milo_obj <- makeNhoods(milo_obj,
                         prop = p, k = k, d = ncomps, reduced_dims = use_rep,
                         refined = TRUE, refinement_scheme = "graph"
  )
  # nhood size
  a <- median(colSums(nhoods(milo_obj)))
  b <- round(mean(colSums(nhoods(milo_obj))), 1)
  d <- range(colSums(nhoods(milo_obj)))
  g <- plotNhoodSizeHist(milo_obj)
  # add text
  plot_data <- ggplot_build(g)$data[[1]]
  x_max <- max(plot_data$x)
  y_max <- max(plot_data$y)
  g <- g +
    annotate("text", label = paste0("mean: ", b), x = x_max * .7, y = y_max * .9, hjust = 0) +
    annotate("text", label = paste0("median: ", a), x = x_max * .7, y = y_max * .8, hjust = 0) +
    annotate("text", label = paste0("range: [", d[1], ", ", d[2], "]"), x = x_max * .7, y = y_max * .7, hjust = 0) +
    labs(title = paste0("k=", k, ", prop=", p), y = "Count") +
    theme_cowplot(
      font_size = 10,
      rel_small = 8 / 10,
      rel_tiny = 8 / 10,
      rel_large = 10 / 10
    )
  return(g)
}


k <- 10
milo.obj <- build_grap(milo.obj, k = k)
p1 <- optimize_k_p(milo.obj, k = k, p = 0.05)
p2 <- optimize_k_p(milo.obj, k = k, p = 0.1)
p <- p1 + p2 +
  plot_layout(ncol = 2)
ggsave2(p,
  filename = file.path(word_dir, DOCNAME, paste0("Fig.k=", k, ".pdf")),
  width = 10, height = 4.5
)

k <- 15
milo.obj <- build_grap(milo.obj, k = k)
p1 <- optimize_k_p(milo.obj, k = k, p = 0.05)
p2 <- optimize_k_p(milo.obj, k = k, p = 0.1)
p <- p1 + p2 +
  plot_layout(ncol = 2)
ggsave2(p,
        filename = file.path(file.path(word_dir, DOCNAME, paste0("Fig.k=", k, ".pdf"))),
        width = 10, height = 4.5
)

k <- 20
milo.obj <- build_grap(milo.obj, k = k)
p1 <- optimize_k_p(milo.obj, k = k, p = 0.05)
p2 <- optimize_k_p(milo.obj, k = k, p = 0.1)
p <- p1 + p2 +
  plot_layout(ncol = 2)
ggsave2(p,
        filename = file.path(file.path(word_dir, DOCNAME, paste0("Fig.k=", k, ".pdf"))),
        width = 10, height = 4.5
)

k <- 30
milo.obj <- build_grap(milo.obj, k = k)
p1 <- optimize_k_p(milo.obj, k = k, p = 0.05)
p2 <- optimize_k_p(milo.obj, k = k, p = 0.1)
p <- p1 + p2 +
  plot_layout(ncol = 2)
ggsave2(p,
        filename = file.path(file.path(word_dir, DOCNAME, paste0("Fig.k=", k, ".pdf"))),
        width = 10, height = 4.5
)


Sys.Date()
sessionInfo()
