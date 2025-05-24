
cross_knn <- function(srat_query = NULL, srat_ref = NULL, k = NULL) {
  # get cells
  coldata.query <- srat_query@meta.data %>%
    rownames_to_column("cell") %>%
    dplyr::select(cell, Origin4, Patient, cellType, cellType2, cellType3)
  coldata.ref <- srat_ref@meta.data %>%
    rownames_to_column("cell") %>%
    dplyr::select(cell, Origin4, Patient, cellType, cellType2, cellType3)

  # get PCA score of each cell
  pc.query <- srat_query@reductions$pca[[]]
  pc.ref <- srat_ref@reductions$pca[[]]

  cat("[INFO] Query KNN ...\n")
  N12 <- queryKNN(pc.ref, query = pc.query, k = k, get.distance = T, BPPARAM = MulticoreParam(workers = 8))
  # Use ascites as reference
  cat("[INFO] Reverse Query KNN ...\n")
  # N21 <- queryKNN(pc.query, query = pc.ref, k = k, get.distance = T, BPPARAM = MulticoreParam(workers = 8))
  cat("[INFO] Find MNN ...\n")
  # mnn.sets <- batchelor::findMutualNN(pc.ref, pc.query, k1 = k, k2 = k, BPPARAM = MulticoreParam(workers = 8))
  N12$cellnames <- rownames(pc.query)
  # N21$cellnames <- rownames(pc.ref)
  knn.list <- list(
    pca.ref = pc.ref, pca.query = pc.query,
    queryToRef = N12, # refToQuery = N21,
    # mnnSets = mnn.sets,
    coldata.ref = coldata.ref, coldata.query = coldata.query
  )
  return(knn.list)
}


construct_crossknn_result <- function(knn.list = NULL) {
  coldata.query <- knn.list$coldata.query
  coldata.ref <- knn.list$coldata.ref
  pc.ref.index <- knn.list$mnnSets$first
  pc.query.index <- knn.list$mnnSets$second
  query.mnn <- unique(pc.query.index)
  query.mnn.ordered <- sort(query.mnn)
  coldata.query[, "nn.cellId1"] <- coldata.ref[knn.list$queryToRef$index[, 1], "cell"]
  coldata.query[, "nn.cellId2"] <- coldata.ref[knn.list$queryToRef$index[, 2], "cell"]
  coldata.query[, "nn.cellId3"] <- coldata.ref[knn.list$queryToRef$index[, 3], "cell"]
  coldata.query[, "nn.cellId4"] <- coldata.ref[knn.list$queryToRef$index[, 4], "cell"]
  coldata.query[, "nn.cellId5"] <- coldata.ref[knn.list$queryToRef$index[, 5], "cell"]
  coldata.query[, "nn.cellId6"] <- coldata.ref[knn.list$queryToRef$index[, 6], "cell"]
  coldata.query[, "nn.celltype"] <- coldata.ref[knn.list$queryToRef$index[, 1], "cellType3"]
  coldata.query[, "nn.tissue"] <- coldata.ref[knn.list$queryToRef$index[, 1], "Origin4"]
  coldata.query[, "nn.donor"] <- coldata.ref[knn.list$queryToRef$index[, 1], "Patient"]
  coldata.query[, "nn.dist"] <- knn.list$queryToRef$distance[, 1]
  coldata.query$cell.anchor <- ifelse(rownames(coldata.query) %in% rownames(knn.list$pca.query)[query.mnn.ordered],
                                      "anchor", "others"
  )
  coldata.query <- coldata.query %>%
    dplyr::mutate(nn.label = paste0(nn.celltype, "@", nn.tissue))
  return(coldata.query)
}

# coldata.mt <- construct_crossknn_result(knn.list = knn.list)
# readr::write_tsv(coldata.mt, file = here::here("output", DOCNAME, t, paste0(cell, "_knn_results.mt2others.tsv")))
################### ==================calculate pearson correlation between each cell pairs from KNN result===============#################


result <- list()
for (i in 1:ncol(data.query)) {
  cell_query <- coldata.query[i, ]$cell
  y <- data.query[, cell_query]
  cell_ref <- queryToRef_idx[i, ]
  ref_mat <- data.ref[, cell_ref]
  cor.mat <- rcorr(as.matrix(y), as.matrix(ref_mat), type = "pearson")
  cor.result <- data.frame(cor.mat$r[1, 2:31])
  colnames(cor.result) <- "cor"
  cor.result$cell <- rownames(cor.result)
  cor.result <- cor.result[order(cor.result$cor, decreasing = T), ]
  pvalue <- data.frame(cor.mat$P[1, 2:31])
  colnames(pvalue) <- "pvalue"
  pvalue$cell <- rownames(pvalue)
  myresult <- plyr::join(cor.result, pvalue, by = "cell")
  result[[i]] <- myresult
}
length(result)




# tissue
t <- "m_pLN"
dir.create(here::here("output", DOCNAME, t), showWarnings = FALSE)









################==========================generate barplot=============##########################
cells.count <- seurat@meta.data %>%
  group_by(Origin4, cellType3) %>%
  dplyr::summarise(cells = n())
cells.count$cellType3 <- as.character(cells.count$cellType3)
df <- inner_join(coldata.mt, cells.count, by = c("nn.tissue" = "Origin4", "nn.celltype_sub" = "cellType3"))
plotdata <- df %>%
  dplyr::mutate(nn.celltype = nn.celltype_sub) %>%
  group_by(nn.tissue, cellType3) %>%
  dplyr::summarise(cells = n()) %>%
  group_by(cellType3) %>%
  dplyr::mutate(frac = cells / sum(cells))

coldata.mt %>%
  dplyr::mutate(nn.celltype = nn.celltype_sub) %>%
  group_by(nn.tissue, cellType3) %>%
  dplyr::summarise(cells = n()) %>%
  group_by(cellType3) %>%
  dplyr::mutate(frac = cells / sum(cells))

p <- ggplot(plotdata, aes(x = nn.tissue, y = frac)) +
  geom_bar(aes(fill = nn.tissue), stat = "identity") +
  geom_text(aes(y = ifelse(frac <= .1, 0.1, frac - .1), label = paste0("N=", cells)), size = 3) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  # scale_fill_manual(values = origin4_color_maps) +
  facet_wrap(. ~ cellType3, ncol = 4) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text.y = element_text(size = 10)
  ) +
  ylab("Fraction") +
  xlab("") +
  ggtitle("Merged Donors")



