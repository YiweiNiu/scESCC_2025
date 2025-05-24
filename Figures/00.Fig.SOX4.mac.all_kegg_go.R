#!/usr/bin/env Rscript

# nohup Rscript analysis/00.Fig.SOX4.mac.all_kegg_go.R > analysis/00.Fig.SOX4.mac.all_kegg_go.log &

DOCNAME <- "00.Fig.SOX4.mac"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

# Tidyverse
library(tidyverse)

# Seurat
library(Seurat)
library(msigdbr)

# anno
library(org.Hs.eg.db)
library(limma)

source(here::here("code/preprocess.R"))

srat = readRDS(here::here('output/04.rm_cells/seurat_mac.rds'))
srat

# KEGG, 359 sets
tab <- getGeneKEGGLinks(species="hsa") %>%
  mutate(PathwayID = str_remove(PathwayID, "path:"))
tab$Symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, tab$GeneID, column="SYMBOL", keytype="ENTREZID")
hsa_pathways <- KEGGREST::keggList("pathway", "hsa") %>%
  tibble(pathway = names(.), description = .) %>%
  mutate(description = str_remove(description, " - Homo sapiens [(]human[])]"))

# get score
gs_lst <- split(tab$Symbol, tab$PathwayID)
srat <- AddModuleScore(srat, features = gs_lst, name = "KEGG", search = TRUE)
gs_score <- srat@meta.data %>%
  rownames_to_column("barcode") %>%
  dplyr::select(barcode, starts_with("KEGG")) %>%
  column_to_rownames("barcode")
colnames(gs_score) <- names(gs_lst)
x <- FetchData(srat, slot = "data", vars = "SOX4")

# cor
y <- cor(x$SOX4, gs_score)
df <- as.data.frame(t(y)) %>%
  rownames_to_column("x") %>%
  dplyr::rename(pathway = x, R = V1) %>%
  left_join(hsa_pathways, by = "pathway")

# top
top50 <- df %>%
  arrange(-R) %>%
  top_n(50, wt = R)
bottom50 <- df %>%
  arrange(R) %>%
  top_n(50, wt = -R)
rbind(top50, bottom50) %>%
  write_csv(file = here::here("output", DOCNAME, "SOX4_cor_top_KEGG.csv"))

# GO bp, 7658 sets
bp_gene_sets = msigdbr::msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
# get score
gs_lst_tmp <- split(bp_gene_sets$gene_symbol, bp_gene_sets$gs_name)
# filt
gs_lst <- list()
#gs_lst_tmp <- gs_lst_tmp[1:10] # test
for(i in 1:length(gs_lst_tmp)) {
  a <- filter_features(srat, gs_lst_tmp[[i]])
  if (length(a) >= 5) { # at least 5 genes
    gs_lst[[names(gs_lst_tmp)[i]]] <- a
  }
}
length(gs_lst)
srat <- AddModuleScore(srat, features = gs_lst, name = "GOBP")
gs_score <- srat@meta.data %>%
  rownames_to_column("barcode") %>%
  dplyr::select(barcode, starts_with("GOBP")) %>%
  column_to_rownames("barcode")
colnames(gs_score) <- names(gs_lst)
x <- FetchData(srat, slot = "data", vars = "SOX4")

# cor
y <- cor(x$SOX4, gs_score)
df <- as.data.frame(t(y)) %>%
  rownames_to_column("x") %>%
  dplyr::rename(pathway = x, R = V1) %>%
  left_join(hsa_pathways, by = "pathway")

# top
top50 <- df %>%
  arrange(-R) %>%
  top_n(50, wt = R)
bottom50 <- df %>%
  arrange(R) %>%
  top_n(50, wt = -R)
rbind(top50, bottom50) %>%
  write_csv(file = here::here("output", DOCNAME, "SOX4_cor_top_GOBP.csv"))


Sys.Date()
sessionInfo()

